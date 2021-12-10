function [y,x,dy_du] = simulateModel(obj,u,varargin)
%SIMULATEMODEL Calculates simulated response y of a NsfPolyModel. 
%The Simulation is based on the Statespace description
% x(k+1) = A*x(k) + B*f(u(k),y(k))
% y(k) = C*x(k) + D*f(u(k))
%where x0 is the initial state of the system. If x0 is not supplied or is empty it will be assumed to be zero.
%The function can further calculate the gradients dy_du which might be helpful for MPC applications.
%If function shall compute Gradients, then the inputnonlinearities f(u,y)
%need to return derivative df/du and df/dy as well. See ./test/test_sim_gradient for examples.
%
% [y,x,dy_du] = simulateModel(obj,u,x0,varargin)
% obj [NsfPolyModel]: The Model to be simulated.
% u [N x nu double]: The inputs columnwise [u(1,:); u(2,:); ... ; u(Ns,:)]
% x0 [nx x 1 double]: The initial state x(1) of the system.
% varargin: Name Value Pair arguments 
%   'CalcInputGradients' [int. vector]: The inputnumbers for which gradients will be computed. 
%                                     	Default: 1:obj.InputDimension 
%   'CalcInputGradientsIxMax' [int scalar]: Maximum time horizon for which gradients will be computed. Def. N 
%   'Ss' [struct]:  Structure containing statespace matrices A, B, C, D to be used for simulation.
%                   If not supplied the these will be extracted from obj.Ss. Default: []
% OUTPUTS:
% y [N x ny matrix]: The outputs columnwise. [y(1,:); y(2,:); ... ; y(Ns,:)]. See addYsample option.
% x [N+1 x nx matrix]: The states columnwise. [x(1,:); x(2,:); ... ; x(Ns+1,:)]
% dy_du [N x ny x nu_ x N_ double]: The gradients dy/du. All dy_du(j,:,:,i) for j<i will be zero since y(j) depends only on u(i) for i<=j. 
%                                   NOTE: nu_ = length(CalcInputGradients), N_ = CalcInputGradientsIxMax  
%                                   1st dim.: Time
%                                   2nd dim.: Outputnumber 
%                                   3rd dim.: Inputnumber corresponding to CalcInputGradients
%                                	4th dim.: Timesample number.
%                                  	dy_du(5,2,3,4) means derivative of dy2(5)/du3(4).

%% Parse ins
% tic;
assert(obj.IsParameterized,'Model is not fully parameterized! Check NoiseVariance and model parameters!')
N = size(u,1);
p = inputParser();
addOptional(p,'X0',[]);
addParameter(p,'Ss',[],@(x) isstruct(x)); 
addParameter(p,'CalcInputGradients',1:obj.InputDimension,@(x) all(x<=obj.InputDimension & mod(x,1) == 0 & x>0 & all(diff(x)>0)));
addParameter(p,'CalcInputGradientsIxMax',N,@(x) x<=N & x>=1);
parse(p,varargin{:});
% a = toc;

%% Init 
% tic;
hasof = obj.HasOutputFeedback;
f = obj.InputNonlinearity;
nv = length(f);
ny = obj.OutputDimension;
if ~isempty(p.Results.Ss)
    Ss = p.Results.Ss; 
    B = Ss.B; D = Ss.D; 
else
    Ss = obj.Ss;
    if nv > 0
        B = Ss.B(:,Ss.InputGroup.Measured); D = Ss.D(:,Ss.InputGroup.Measured); 
    else
        B = zeros(size(Ss.A,1),nv); D = zeros(ny,nv); 
    end
end
A = Ss.A; C = Ss.C;   nx = size(A,1); 
alpha = cell(1,nv); % get Parametervectors of Input Nonlinearity
for ni = 1:nv; alpha{ni} = obj.getFiParams(ni); end
f_nouts = obj.F_nargout;                    % Number of output arguuments of Input Nonlinearity
ix_no_of = find(~hasof);                    % Inputnonlin without outputfeedback
ix_of = find(hasof);                        % Inputnonlin with outputfeedback
ix_lin = arrayfun(@(fi) isempty(fi.fun),f);	% linear input f(u,y) = u
u_grads = p.Results.CalcInputGradients;
Ndu = p.Results.CalcInputGradientsIxMax;
x0 = p.Results.X0(:);
if isempty(x0); x0 = zeros(nx,1); end
% b = toc;

%% Malloc
% tic;
y = NaN(ny,N); x = NaN(nx,N+1); x(:,1) = x0;
fuy = NaN(nv,1);
if nargout>2
    dx_du = zeros(nx,length(u_grads),N,Ndu); %dx(k)/du(i)
    dy_du = zeros(ny,length(u_grads),N,Ndu);
    df_du = zeros(nv,length(u_grads));
    df_dy = zeros(nv,ny);
    ix_du = false(nv,length(u_grads)); %index vector
    ix_grad = cell(nv,1);
    need_grad = false(nv,1);
    for j = 1:nv
        if isempty(f(j).fun)  %assign 1 to f_j = u_input_idx
            df_du(j,f(j).input_idx==u_grads) = 1;
        else
            ix_du(j,:) = arrayfun(@(ix_u) any(f(j).input_idx == ix_u),u_grads);
        end
        %ix_grad{j} = arrayfun(@(ix) any(ix == u_grads),f(j).input_idx);
        ix_grad{j} = cell2mat(arrayfun(@(ix) find(f(j).input_idx==ix),u_grads,'UniformOutput',0));
        need_grad(j) = any(arrayfun(@(ix) any(ix == u_grads),f(j).input_idx));
    end
end

%% Simulate: x(k+1) = Ax(k) + B*f(u(k),y(k))
for k = 1:N
    % Evaluate Inputnonlinearities without outputfeedback
    for ni = ix_no_of'
        if ix_lin(ni) % linear in
            fuy(ni,1) = u(k,f(ni).input_idx);
        else
            out = cell(1,f_nouts(ni)-1);
            [fuy(ni,1), out{:}] = f(ni).fun(u(k,f(ni).input_idx),alpha{ni});
            if nargout > 2 && need_grad(ni)
                df_du(ni,ix_du(ni,:)) = out{2}(ix_grad{ni});
            end
        end
    end
    
    % First Compute output @k
    y(:,k) = C*x(:,k) + D(:,~hasof)*fuy(~hasof,1) + obj.OutputOffset; %[y(1) y(2) ... ]

    % Evaluate Inputnonlinearities with outputfeedback
    for ni = ix_of' 
        out = cell(1,f_nouts(ni)-1);
        [fuy(ni,1), out{:}] = f(ni).fun(u(k,f(ni).input_idx),alpha{ni},y(f(ni).output_idx,k)');
      	if nargout > 2 && need_grad(ni)
            df_du(ni,ix_du(ni,:)) = out{2}(ix_grad{ni});
            df_dy(ni,f(ni).output_idx) = out{3};
        end
    end
    % Now compute new state of the system
    x(:,k+1) = A*x(:,k) + B*fuy; %[x(1) x(2)]

    if nargout>2 % Compute Gradient if necessary 
        for i=1:min(k,Ndu)
            % Compute dy(k)/du(i)
            dy_du(:,:,k,i) = C*dx_du(:,:,k,i); % []
            if i == k
                dy_du(:,:,k,i) = dy_du(:,:,k,i) + D(:,~hasof)*df_du(~hasof,:);
            end

            % Comp. dx(k+1)/du(i) = A dx(k)/du(i) + B df(u(k),y(k))/du(i)
            dx_du(:,:,k+1,i) = A*dx_du(:,:,k,i) + B*df_dy*dy_du(:,:,k,i);
            if i==k 
                dx_du(:,:,k+1,i) = dx_du(:,:,k+1,i) + B*df_du; 
            end
        end
    end
end
x = x'; 
y = y';
if nargout>2
    dy_du = permute(dy_du,[3 1 2 4]);
end
% c = toc;
% fprintf('a=%f ; b=%f; c=%f \n',a/(a+b+c),b/(a+b+c),c/(a+b+c));
end

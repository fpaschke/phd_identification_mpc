function yp = calcPredictions(obj,y,varargin)
%Calculate Multi-Step-Ahead predictions y[t+k|t] (k = 1,2,...,Hp) of NsfPolyModel at each time sample t. 
%The simulation is similar to simulateModel and is based on the statespace description
% x[t+1] = A*x[t] + B*f(u[t],y[t])
% y[t] = C*x[t] + D*f(u[t])
%The states x[t] are estimated within the function by running a kalman-filter on the measured sequence u and y. 
%This function is a vectorized implementation of simulateModel, such that the simulation can be done
%for multiple initial states x[t]. Thus it is much faster then calling simulate for each state x[t] seperately.
%This method is meant to be used for model validation. 
%
% yp = calcPredictions(obj,y,u,varargin)
% obj [NsfPolyModel]: An identified model.
% y [(cell of) N x ny double]:  The Measured Outputs of the Model. Each column represents one output.
% u [(cell of) N x nu double]:  The Measured Inputs. Each column represents one input.
% yp [Nset cell of N x ny x Hp double]:	Each column of y represents a prediction.  yp{1}(:,1,1) represents the 
%                                   	1-step-ahead predictions of output 1 (dataset 1), where yp{3}(:,2,10) represents 
%                                    	the 10-step-ahead predictions of output 2 (dataset 3).
% varargin: Name Value Pairs
%   'Hp' [opt. pos int]:   	Prediction Horizon in #Samples (default: 1)
%   'Nmax' [opt. pos int]:	Maximum #Samples from which the kalman filter will be initialized (default: 'auto').
%   'UseGls' [opt. log.]:	If true then the initial conditions will be calculated by GLS. 
%                           GLS is theoratically more accurate then oridinary LS. 
%                           However the difference in most cases is not significant. Def. false


%% PARSE INS
% handle old interface: yp = calcPredictions(obj,y,u,Hp,Hp_obs,Init_oe,x0obs)
if length(varargin)>1 && ~ischar(varargin{2}) 
    u = varargin{1};
    if isempty(varargin{2}); Hp = 1; else; Hp = varargin{2}; end
    if length(varargin)<3 || isempty(varargin{3}); Hp_obs = 'auto'; else; Hp_obs = varargin{3}; end
    if length(varargin)<4 || isempty(varargin{4}); UseGls = true; else; UseGls = varargin{4}; end
    if length(varargin)<5 || isempty(varargin{5}); X0 = []; else; X0 = varargin{5}; end
    wrn = true;
else
    p = inputParser(); p.KeepUnmatched = true;
    addOptional(p,'u',[]);
    addParameter(p,'Hp',1,@(x) x>0 && isscalar(x));
    addParameter(p,'Nmax','auto');
    addParameter(p,'UseGls',false,@(x) x == 1 || x == 0);
    addParameter(p,'Warnings',true,@(x) x == 1 || x == 0);
    addParameter(p,'X0',[]);
    parse(p,varargin{:});
    Hp = p.Results.Hp; Hp_obs = p.Results.Nmax; UseGls = p.Results.UseGls; X0 = p.Results.X0; u = p.Results.u; wrn = p.Results.Warnings;
end
   
%% Treat Wrn
if ~wrn
    wrnStruct = warning;
    warning('off');
end

%% PREPARE DATA AND INIT
sys = obj.Ss;
assert(~isempty(sys),'System isnt fully parameterized which means that any of the polynomials A,B,... or E contains NaNs!');
if ~iscell(y); cellin = false; else; cellin = true; end
[y,u] = obj.getRawData(y,u);
Nsets = length(y);

%% SIMULATE
assert(Hp>0 && mod(Hp,1)==0,'Hp needs to be a pos. integer!');
if any(obj.HasOutputFeedback); warning('Using Deterministic Model for prediction!'); end
if isempty(X0)
    X0 = obj.calcIc(y,u,'Nmax',Hp_obs,'UseGls',UseGls); % Initial values for observer
end
fu = obj.evalInputNonlinearity(u,y); % Evaluate InputNonlinearity
[Yobs, Xobs, Eobs] = idModels.alg.lti.observe(sys,y,fu,X0,obj.OutputOffset); %xobs(1|0) = X0   ... %yobs(1)=ypred(1) ...

yp = cell(length(y),1);
for ns = 1:Nsets 
    yp{ns,1} = simnl(obj,sys,u{ns,1},Xobs{ns,1},Hp);
end
if ~cellin; yp = yp{1}; end

%% Restore Warning state
if ~wrn
    warning(wrnStruct);
end
end

%% vectorized simulation for nonlinear system
function y = simnl(obj,sys,u,xo,Hp) 
    % Init
  	A = sys.A; C = sys.C;  nu = size(obj.B,2);
    if nu > 0
        B = sys.B(:,sys.InputGroup.Measured); D = sys.D(:,sys.InputGroup.Measured);
    else
        B = []; D = [];
    end
    N = size(xo,1);
    ny = size(C,1);
    y = NaN(N,ny,Hp);
    y0 = obj.OutputOffset;
    f = obj.InputNonlinearity;
    nv = length(f) ;
    % Compute Inputnonlinearities without feedback
    k = 1; ix = find(~obj.HasOutputFeedback);
    fu = NaN(N,length(ix)); %inputs to linear block of the system without outputfeedback
   	
    hasof = obj.HasOutputFeedback;
    f_nouts = obj.F_nargout; % Number of output arguuments of Input Nonlinearity
    ix_no_of = find(~hasof); ix_of = find(hasof);
    fuy = NaN(N-k+1,nu);
	alpha = cell(1,nv); % get Parametervectors of Input Nonlinearity
    for ni = 1:nv
        alpha{ni} = obj.getFiParams(ni);
    end
    for ni = ix_no_of' % "pre-evaluate" all inuputnonlinearities without outputfeedback
        out = cell(1,f_nouts(ni)-1); 
        if ~isempty(f(ni).fun)
            [fuy(:,ni), out{:}] = f(ni).fun(u(1:N-k+1,f(ni).input_idx),alpha{ni});
        else
            fuy(:,ni) = u(1:N-k+1,f(ni).input_idx);
        end
    end
    
    % Do Simulation
    x = xo'; %states at current simulation step % [x(1|0) x(2|1) ... x(N|N-1)]
    for k = 1:Hp
        % Compute output of k-th time step (k-step ahead prediction)
        yy  = C*x; %k:end oder 
        if nu > 0 
            yy = yy + D(:,~hasof)*fuy(k:N,~hasof)'; %dim1: ny dim2: t
        end
        % Assign Output
        y(k:end,:,k) = bsxfun(@plus,yy',y0'); %dim1: t %dim2: y %dim3: k/Hap
        % Compute InputNonlinearities
        for ni = ix_of'
            out = cell(1,f_nouts(ni)-1);
            [fuy(k:end-1,ni), out{:}]  = f(ni).fun(u(k:end-1,f(ni).input_idx),alpha{ni},y(k:end-1,f(ni).output_idx,k));
        end
        % Now compute the new states of the system
        x = A*x(:,1:end-1); % [x(2|0) x(3|1) ... x(N+1|N-1)]
        if nu > 0 
            x = x + B*fuy(k:end-1,:)'; 
        end
    end
end
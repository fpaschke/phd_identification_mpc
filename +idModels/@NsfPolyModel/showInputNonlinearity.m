function fig = showInputNonlinearity(obj,ix_f,u0,varargin)
%SHOWINPUTNONLINEARITY Plots  
%
% f = showInputNonlinearity(obj,ix_f,u0,y0,varargin)
%
% obj [NsfPolyModel]:	Instance of NsfPolyModel
% ix_f [pos. int]:      Index of component of f that should be shown, i. e. f(ix_f) (will be plotted to y-axis). 
% u0 [1 x nv vector]:   nv = length(obj.InputNonlinearity(ix_f).input_idx)  	
%                       The point u0 that will be evaluated. One or two elements can be NaN. 
%                       The NaN input(s) will be varied from min to max (obj.InputMin and obj.InputMax)!
%                    	The values that are not NaN specify the remaining inputvalues!  
%                      	
% y0 [1 x nv vector]:	nv = length(obj.InputNonlinearity(ix_f).output_idx)  	
%                       The point y0 that will be evaluated. One or two elements can be NaN. 
%                       The NaN outputs(s) will be varied from min to max (obj.InputMin and obj.InputMax)!
%                    	The values that are not NaN specify the remaining inputvalues!
% varargin: Optional name Value pairs:
% 'u' [N x 1 double]:   Points at which f(ix_in) will be evaluated. Def.: linspace(InputMin(ix_in),InputMax(ix_in),1e3)
% 'Axes' [Axes]:        If specified then Inputnonlinearity will be plotted to specific axes object
% 'Confidence' [opt. pos scalar double 0...1]:  Confidence bound (Def. 0.99). .99 means that 99% of estimates will be within
%                                               that bound. Higher values lead to bigger bounds. 1 leads to Inf.
%
% EXAMPLES:
% rng(0); clear; u = [kron(rand(1,20),ones(1,50))' kron(rand(1,50),ones(1,20))' 1+kron(rand(1,10),ones(1,100))']; e = .1*randn(size(u,1),1);
% F(1).fun = @idModels.func.fun_powdT; F(1).parameters = [2 1]; F(1).free = [1 0]; F(1).input_idx = [1 3]; F(1).output_idx = 1;
% F(2).fun = @idModels.func.fun_logistic; F(2).parameters = [.5 10]; F(2).free = [1 0]; F(2).input_idx = 2;
% A = conv([1 -.9],[1 -.95]); C = conv([1 -.7+.4i],[1 -.7-.4i]); B1 = [0 .1]; B2 = [0 .2]; 
% M = idModels.NsfPolyModel(2,[1 1],[1 1],2,'InputNonlinearity',F);
% y(1:2,1) = 0; 
% for k = 3:length(u)
%   [f1,~] = F(1).fun(u(k-1,F(1).input_idx),F(1).parameters,y(k-1,F(1).output_idx)); 
%   [f2,~] = F(2).fun(u(k-1,F(2).input_idx),F(2).parameters); 
%   y(k,1) = -A(2:end)*y(k-1:-1:k-2) + B1(2:end)*f1 + B2(2:end)*f2 + C*e(k:-1:k-2); 
% end
% M.identify(y,u);
% M.show('Prediction',y,u,'Hp',10)
% figure('color','w');
% draw f1(u1,u3=1) and f1(u1,u3=2) @ y=0
% s1 = subplot(1,3,1); 
% M.showInputNonlinearity(1,[NaN 1],0,'Axes',s1,'Confidence',.99,'Color','b'); hold on;
% M.showInputNonlinearity(1,[NaN 2],0,'Axes',s1,'Confidence',.99,'Color','r'); 
% % draw f(u1,u3) @ y=0
% s2 = subplot(1,3,2); M.showInputNonlinearity(1,[NaN 0],NaN,'Axes',s2,'Confidence',.99);
% % draw f2(u=[u1]) @ y=0
% s3 = subplot(1,3,3); M.showInputNonlinearity(2,[NaN],'Axes',s3,'Confidence',.99);%% CHECK ix_f,u0

assert(obj.HasInputNonlinearity(ix_f),['The obj.InputNonlinearity(' num2str(ix_f) ') is identity!']);
assert(isscalar(ix_f) && ix_f<=length(obj.InputNonlinearity),'ix_f needs to be <= length(obj.InputNonlinearity)');
assert((isvector(u0) && length(u0) == length(obj.InputNonlinearity(ix_f).input_idx) && sum(isnan(u0))==1) || ...
       (ismatrix(u0) && size(u0,2) == length(obj.InputNonlinearity(ix_f).input_idx)),...
            'u0 needs to be vector of length length(obj.InputNonlinearity(ix_f).input_idx) of which one element is NaN or a matrix with as many columns as length(obj.InputNonlinearity(ix_f).input_idx)!');
    
%% PARSE varargin
p = inputParser();
addOptional(p,'y0',zeros(obj.InputNonlinearity(ix_f).output_idx,1),@(x) length(x) == length(obj.InputNonlinearity(ix_f).output_idx));
addParameter(p,'Confidence',.99,@(x) x>=0 && x<=1);
addParameter(p,'Color',lines(1));
addParameter(p,'Axes',[]);
parse(p,varargin{:});

%% INIT
covPk = obj.Info.CovP;
if isnumeric(covPk) && ~any(isnan(covPk(:)))
    plotConf = true;
else
    plotConf = false;
	warning('Couldnt compute confidence bounds because parametercovariance (obj.Info.CovP) is not available!');
end

% Find varying ins 
N = 25; 
if isvector(u0)
    ix_in = find(isnan(u0));
    u0 = u0(:)';
else
    error('Invalid input for u0!');
end

% Find varing out
y0 = p.Results.y0;
if isvector(y0)
    ix_out = find(isnan(y0));
    y0 = y0(:)';
else
    error('Invalid input for y0!');
end
assert(~isempty([ix_in ix_out]) > 0 && length([ix_in ix_out]) < 3,'You need to specify the components of u (and y) that should be varied using NaN! Seed documentation!. The Maximum number of NaNs (varying parameters) are 2!');

%% Build u and y Evaluate f(u,p,...)
% Build u and y
u = u0.*ones(N,length(obj.InputNonlinearity(ix_f).input_idx));
u(:,ix_in) = cell2mat(arrayfun(@(i) linspace(obj.InputMin(i),obj.InputMax(i),N)',obj.InputNonlinearity(ix_f).input_idx(ix_in),'UniformOutput',0));
y = y0.*ones(N,length(obj.InputNonlinearity(ix_f).output_idx));
y(:,ix_out) = cell2mat(arrayfun(@(i) linspace(obj.OutputMin(i),obj.OutputMax(i),N)',obj.InputNonlinearity(ix_f).output_idx(ix_out),'UniformOutput',0));

if length([ix_in ix_out]) == 1; NN = 1; else; NN = N; end
f = NaN(N,NN);
confF = f;
out = cell(1,obj.F_nargout(ix_f)-2);
U = u; Y = y;
for n = 1:NN
    if length(ix_in) == 2 % both varying nondependent vars are u's
        U(:,ix_in(2)) = u(n,ix_in(2)); 
    elseif length(ix_out) == 2 % both varying nondependent vars are y's
        Y(:,ix_out(2)) = y(n,ix_out(2)); 
    elseif length(ix_in) == 1 && length(ix_out) == 1 % varying nondependent are u and y
        Y(:,ix_out(1)) = y(n,ix_out(1)); 
    end
    if obj.HasOutputFeedback(ix_f)
        [f(:,n),df,out{:}] = obj.InputNonlinearity(ix_f).fun(U,obj.getFiParams(ix_f),Y);
    else
        [f(:,n),df,out{:}] = obj.InputNonlinearity(ix_f).fun(U,obj.getFiParams(ix_f));
    end

    %% Calc confidence Bound
    if plotConf
        if size(df,2) ~= sum(obj.InputNonlinearity(ix_f).free) % returns df/dp wrt to full parameter vector
            df = df(:,obj.InputNonlinearity(ix_f).free);
        end
        Npar = size(covPk,1);
        dF = zeros(N,Npar);

        % find corresponding index in Pvec
        ix = Npar - sum(arrayfun(@(fi) sum(fi.free),obj.InputNonlinearity(end:-1:ix_f))) + 1;
        dF(:,ix:ix+sum(obj.InputNonlinearity(ix_f).free)-1) = df;
        stdF = sqrt(diag(dF*covPk*dF'));
        confF(:,n) = util.norminv((p.Results.Confidence+1)/2,zeros(N,1),stdF);
    end
end

%% DRAW
if isempty(p.Results.Axes); figure('color','w'); end
if length([ix_in ix_out]) == 1 % 2D
    if ~isempty(ix_in)
        x = u(:,ix_in); xlab = obj.getInLab(obj.InputNonlinearity(ix_f).input_idx(ix_in));
    else
        x = y(:,ix_out); xlab = obj.getOutLab(obj.InputNonlinearity(ix_f).output_idx(ix_out));
    end
    plot(x,f,'LineStyle','-','Color',p.Results.Color); hold on;
  	if plotConf; plot(x,f+[confF -confF],'LineStyle','--','Color',p.Results.Color); end
    ylab = ['Eingangsnichtlinearit\"at $f_{' num2str(ix_f) '}$'];
    zlab = [];
else
    if length(ix_in) >= 1
        y_ = u(:,ix_in(1)); ylab = obj.getInLab(obj.InputNonlinearity(ix_f).input_idx(ix_in(1)));
    else
        y_ = y(:,ix_out(1)); ylab = obj.getOutLab(obj.InputNonlinearity(ix_f).output_idx(ix_out(1)));
    end
  	if length(ix_in) == 2 % both varying nondependent vars are u's
        x = u(:,ix_in(2)); xlab = obj.getInLab(obj.InputNonlinearity(ix_f).input_idx(ix_in(2)));
    elseif length(ix_out) == 2 % both varying nondependent vars are y's
        x = y(:,ix_out(2)); xlab = obj.getOutLab(obj.InputNonlinearity(ix_f).output_idx(ix_out(2)));
    elseif length(ix_in) == 1 && length(ix_out) == 1 % varying nondependent are u and y
        x = y(:,ix_out(1)); xlab = obj.getOutLab(obj.InputNonlinearity(ix_f).output_idx(ix_out(1)));
    end
    
    mesh(x,y_,f); hold on
    if plotConf
        mesh(x,y_,f+confF,'LineStyle',':'); hold on
        mesh(x,y_,f-confF,'LineStyle',':');
    end
    zlab = ['Eingangsnichtlinearit\"at $f_{' num2str(ix_f) '}$'];
end
util.formatFigure(12,xlab,ylab,zlab,[],[],'Axes',gca); 
fig = gcf;
end
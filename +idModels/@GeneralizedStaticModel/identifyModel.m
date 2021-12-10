function e = identifyModel(obj,y,u,opt)
%Performs identification of GeneralizedStaticModel.
%
% e = identifyModel(obj,y,u,opt)
% obj [PolynomialModel]: Any subclass of BModel with an identifyModel method.
% y [cell of N x ny double Matrix]: Raw Output Data for Model identification
% u [cell of N x nu double Matrix]: Raw Input Data for Model identification
% e [cell of N x nu double Matrix]: Residuals of identifed Model.
% opt [struct]: Structure Specifying the options to use for estimation -> see identifyModelOptions

%% IMPORT AND CHECKS
import idModels.func.*;
assert(obj.OutputDimension == 1,'Only Single-Output Models supported yet!');
ls_solvers = {'lsqnonlin' 'LEVMAR' 'LM_DER' 'MKLTRNLS' 'NL2SOL'};

%% REMOVE OUTPUT OFFSET
if ~all(obj.OutputOffset==0)
    y = cellfun(@(yi) yi - obj.OutputOffset',y,'UniformOutput',0);
end

%% INIT CELL2MAT
Ns = cellfun(@length,y);
y = cell2mat(y);
u = cell2mat(u);
N = size(y,1);

%% GET Startingvalue of Parameter
[p0,~,p] = obj.getPvec;
nP = length(p0);
if opt.EstimateOutputOffset
   p0 = [p0(:); 0]; % Da OutputOffset in identify subtrahiert wird
end
assert(isnumeric(p0) && all(~isnan(p)) && ~isempty(p0),'You need to supply a initial guess for the Parameters of the Model');

%% if deal is used 
for no = 1:obj.OutputDimension
    f{no,1} = str2func(obj.Fun{no});
    hasDeal(no,1) = ~isempty(strfind(obj.Fun{no},'deal'));
end

%% Do Optimization
Unmatched = util.struct2namevaluepairs(opt.Unmatched);
ix = regexp(opt.Solver,'_');
tb = opt.Solver(1:ix-1);
solver = opt.Solver(ix+1:end);
if strcmpi(opt.Cost,'abs') && strcmpi(solver,'lsqnonlin')
   warning('You cannot use a Least squares solver for absolute value problems! fminunc without gradients will be used instead!')
   tb = 'matlab';
   solver = 'fminunc';
   opt.UseGradient = false;
end

if strcmpi(tb,'matlab')
	options = optimoptions(solver,'Display',opt.Display,'SpecifyObjectiveGradient',opt.UseGradient,...
                        'CheckGradients',opt.CheckGradient,Unmatched{:});
    if strcmpi(opt.Cost,'quad')
        popt = lsqnonlin(@resid,p0,[],[],options);
    else
        popt = fminunc(@absValFit,p0,options);
    end
elseif strcmpi(tb,'opti')
    if opt.CheckGradient opt.CheckGradient = 'on'; else opt.CheckGradient = 'off'; end
    options = optiset('solver',lower(solver),'display',opt.Display,'derivCheck',opt.CheckGradient,Unmatched{:});
    args = {'x0' p0 'options' options};
    if any(strcmpi(solver,solvers)) % since NLS Solvers expect y and y_hat we need to compute ydata    
        args = [{'fun' @predict 'ydata' y(:)/sqrt(N)} args];
        if opt.UseGradient args = [args {'grad' @opti_dy_grad}]; end
    else
        error('Unkonown Solver')
    end
	Opt = opti(args{:});
    popt = solve(Opt);
else
    error('Unkonwon Solver/Toolbox')
end

%% Calc residuals
e_mat = resid(popt)*sqrt(N); % residuals
e = mat2cell(e_mat,Ns,ones(1,obj.OutputDimension)); % Put to cell format again
obj.NoiseVariance = obj.calcCovariance(e,length(popt));

%% Set optimal parameters and write Info struct
obj.setPvec(popt(1:end-sum(opt.EstimateOutputOffset)));
y0 = popt(nP+1:end); 
if any(opt.EstimateOutputOffset)
    obj.OutputOffset(logical(opt.EstimateOutputOffset)) = obj.OutputOffset(logical(opt.EstimateOutputOffset)) + y0;
end

Info = obj.Info; 
Info.NumEstParameters = length(popt);   
obj.Info = Info;

%% NESTED FUNCTIONS
function J = absValFit(P)
  	e = resid(P);
    J = sum(abs(e/sqrt(N)));
end
    
function [e, de] = resid(P)
    if nargout > 1
        [yp, dy] = predict(P);
        de = -dy; 
    else
        yp = predict(P);
    end
    e = y/sqrt(N) - yp;
    e = e(:);
end

function [yp, dy] = predict(P)
    obj.setPvec(P(1:end-sum(opt.EstimateOutputOffset)));
    old_p = P;
    for no = 1:obj.OutputDimension
        if nargout > 1 || hasDeal(no) 
            [yp(:,no),dy{no}] = f{no}(u,obj.Parameters{no});
            if size(dy{no},2) ~= sum(obj.Free{no}) % Function seems to return gradient wrt to full parameter vector!
                dy{no} = dy{no}(:,obj.Free{no});
            end
            dy0{no} = ones(N,opt.EstimateOutputOffset(no));
        else
            yp(:,no) = f{no}(u,obj.Parameters{no});
        end
        if opt.EstimateOutputOffset(no)
            yp(:,no) = yp(:,no) + P(nP+sum(opt.EstimateOutputOffset(1:no)));
        end
    end
    yp = yp(:)/sqrt(N);
    if nargout > 1
        dy = [blkdiag(dy{:}) blkdiag(dy0{:})]/sqrt(N);
        old_dy = dy;
    end
end

% Gradient for optis lsqnonlin
function dyp = opti_dy_grad(P)
    if any(P ~= old_p) || isempty(old_dy)  % we need to recompute the gradient :-/
        [~, dyp] = predict(P);
    else % we can save some time here :-)
        dyp = old_dy;
    end
end


end
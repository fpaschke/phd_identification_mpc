function e = identifyModel(obj,y,u,opt)
%Performs identification of PolynomialModel.
%
% e = identifyModel(obj,y,u,opt)
% obj [PolynomialModel]: Any subclass of BModel with an identifyModel method.
% y [cell of N x ny double Matrix]: Raw Output Data for Model identification
% u [cell of N x nu double Matrix]: Raw Input Data for Model identification
% e [cell of N x nu double Matrix]: Residuals of identifed Model.
% opt [struct]: Structure Specifying the options to use for estimation -> see identifyModelOptions

assert(obj.OutputDimension == 1,'Identification for MultiOutput Systems is unsupported yet!');

%% REMOVE OUTPUT OFFSET
if ~all(obj.OutputOffset==0)
    y = cellfun(@(yi) yi - obj.OutputOffset',y,'UniformOutput',0);
end

%% CELL2MAT
Ns = cellfun(@length,y);
N = sum(Ns);
y = cell2mat(y);
u = cell2mat(u);

%% NORMALIZE INPUTS
u_norm = max([abs(obj.InputMin)'; abs(obj.InputMax)']);
u = bsxfun(@times,u,1./u_norm);
y_norm = max([abs(obj.OutputMin)'; abs(obj.OutputMax)']);
y = bsxfun(@times,y,1./y_norm);

%% Remove fixed part from y
[~,ixfree,allP] = obj.getPvec(u_norm,y_norm);
PHI = cell2mat(obj.getRegressorMatrices(u));
if any(~ixfree)
    y = y - PHI(:,~ixfree)*allP(~ixfree);
end
assert((opt.FirstDerivative==0 && opt.SecondDerivative==0) || all(ixfree),'Fixed parameters and constraints cant be handled now!');


%% GET REGRESSOR MATRIX OF FREE PARAMS
PHI = PHI(:,ixfree);
if opt.EstimateOutputOffset
   PHI = [ones(N,1) PHI]; 
end

%% GET REGRESSOR MATRICES FOR CONSTRAINTS
if ~(opt.FirstDerivative==0) || ~(opt.SecondDerivative==0)
    uv = cell2mat(arrayfun(@(ni) linspace(obj.InputMin(ni)./u_norm(ni),obj.InputMax(ni)./u_norm(ni),opt.Nconstraints)',1:obj.InputDimension,'UniformOutput',0));
    [~,dPhi,ddPhi] = obj.getRegressorMatrices(uv);
end

%% ESTIMATE
if opt.FirstDerivative==0 && opt.SecondDerivative==0
    P = PHI\y; %estimate
else
    A = []; b = [];
    if opt.FirstDerivative==1 % Constraints regarding 1st derivative 
        A = [A; -blkdiag(dPhi{:})]; b = [b; zeros(obj.InputDimension*opt.Nconstraints,1)];
    elseif opt.FirstDerivative==-1
        A = [A; blkdiag(dPhi{:})]; b = [b; zeros(obj.InputDimension*opt.Nconstraints,1)];
    end
   	if opt.SecondDerivative==1 % Constraints regarding 2nd derivative  
        A = [A; -blkdiag(ddPhi{:})]; b = [b; zeros(obj.InputDimension*opt.Nconstraints,1)];
    elseif opt.FirstDerivative==-1
        A = [A; blkdiag(ddPhi{:})]; b = [b; zeros(obj.InputDimension*opt.Nconstraints,1)];
    end
    A = [zeros(size(A,1),opt.EstimateOutputOffset) A];
    P = quadprog(2*(PHI'*PHI),-2*y'*PHI,A,b);
end

%% CALC RESIDUALS
e_mat = y_norm*(y - PHI*P); % residuals
e = mat2cell(e_mat,Ns,ones(1,obj.OutputDimension)); % Put to cell format again
obj.NoiseVariance = obj.calcCovariance(e,length(P));

%% WRITE TO OBJ
Info = obj.Info; 
Info.NumEstParameters = length(P); 
Info.covP = e_mat'*e_mat/N.*inv(PHI'*PHI);    
obj.Info = Info;
y0 = y_norm*P(1);
if any(opt.EstimateOutputOffset)
    obj.OutputOffset(logical(opt.EstimateOutputOffset)) = obj.OutputOffset(logical(opt.EstimateOutputOffset)) + y0;
    P(1) = []; 
end
obj.setPvec(P,u_norm,y_norm);
end
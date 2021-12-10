function [Cost,M] = estModel(y,u,varargin)
% Instantiates an NsfPolyModel with given Parameters and performs Identification and Testing. 
% This function can be used for model order selection.
%
% [Cost,M] = estFun(y,u,varargin)
% y [(cell of) N x ny double Matrix]: Output Data for Model identification
% u [(cell of) N x nu double Matrix]: Input Data for Model identification
% Cost [1 x 3 double]: RMSE or MAE (if 'UseMae'= true) of [HpVal-Multistep HpVal-Step 1-Step] prediction errors.
% M [NsfPolyModel]: The identified Model.
% varargin: Opt. Name Value Pairs:
% MODEL SETUP AND IDENTIFICATION:
%   'Type' [char]:                      Parametrization of the model. Def.: 'ARMAX'
%                                       Supported: 'ARX' 'BJ' 'ARMAX' 'OE'. Def. 'ARMAX'. 
%   'Order' [pos. int scalar]:          Order of A,B,C,D,E- Polynomials. Def. 1
%   'SeasonalityC' [pos. int scalar]:	Adds a free ciefficient c_s to C polynomial, i. e.
%                                       C(q^-1) = Cn(q^-1) + c_s q^-s (Cn(q) is C polynomial of ARMAX/BJ structure)
%                                       (Cn(q-^1) = 1 + c_1 q^-1 + ... c_n q^-n). Def.: []
%   'SeasonalityD' [pos. int scalar]:	Adds a free ciefficient d_s to D polynomial, i. e.
%                                       D(q^-1) = Dn(q^-1) + d_s q^-s (Dn(q) is D polynomial of BJ structure)
%                                       (Dn(q-^1) = 1 + d_1 q^-1 + ... d_n q^-n). Def.: []
%   'SeasonalityOrder' [int scalar>=1 or 'ByOrder']:       
%                                       If 'ByOrder' then seasonal part of C and D will be increased with parameter Order n, i. e. 
%                                       D(q-^1) = Dn(q-^1) +  d_(s-n+1)q^(-(s-n+1)) + ... + d_sq^-s + ... + d_(s+n-1)q^(-s+n+1)
%                                       (C(q^-1) analog).
%                                       If integer then n in D(q-^1) will be a fixed value specified by this parameter
%   'MdlProp' [struct]:                	Struct of Properties to be set for the Model structure. Def. struct()
%                                       E. g. MdlProp.Ts = 5 and MdlProp.TimeUnit = 'min' sets sampling 
%                                       time to 5 and TimeUnit to 'min'
%   'IdOpt' [struct]:                   Struct of Options to be passed to identify(y,u,...). Def.: struct()
%                                       IdOpt.Hp = 10 set identification horizon to 10.
%   'InitFromPreviousN' [logical]:      If true function will save the previously estimated model in an persistent variable that
%                                       will be used for initialization in next estimation if possible. This works if the 
%                                       previously estimated model is contained in the model to be estimated. If not the standard
%                                       initialization method IdOpt.InitMethod will be used. 
%                                       This option is meant to be used in combination with idModels.test.sweepPar(...), to try 
%                                       different model orders, then make sure that first parameter to be sweeped is 'Order', e. g. 
%                                       idModels.test.sweepPar({'Order' 'Type'},@idModels.NsfPolyModel.estModel,y,u,'InitFromPreviousN',true,'Order',1:10,'Type',{'arx' 'armax' 'bj' 'oe'})
%   'FactorizeG' [logical]:             If true the factorized parametrization of A,B and E will be used. Def.: false
%   'MinRealTol' [double]:           	Performs Model order reduction of Statespace object using Matlabs minreal(ss,MinRealTol). Def.: [] (no order reduction)
% VALIDATION:
%   'HpVal' [pos. int scalar]:          Prediction Horizon used for validation.
%   'UseMae' [logical]:                 Outputs MMAE instead of RMMSE.
%
% rng(0); u = kron(randn(20,1),ones(50,1)); e = .1*randn(length(u),1);
% A = conv([1 -.9],[1 -.8]); B = [0 1]; C = [1 .8];
% y = filter(B,A,u) + filter(C,A,e);
% res1 = idModels.NsfPolyModel.estModel(y,u,'Order',1,'Type','ARMAX','HpVal',10);
% res2 = idModels.NsfPolyModel.estModel(y,u,'Order',2,'Type','ARMAX','HpVal',10);
% res3 = idModels.NsfPolyModel.estModel(y,u,'Order',3,'Type','ARMAX','HpVal',10);
% bar([res1(1) res2(1) res3(1)]); xlabel('Order of Polynomials n'); ylabel('MMSE of 10-Step-Prediction-Errors')

%% PARSE INS
p = inputParser; p.KeepUnmatched = true;

% MODEL AND IDENT
addParameter(p,'Type','ARMAX',@(x) any(strcmpi(x,{'ARX' 'BJ' 'ARMAX' 'OE'}))); 
addParameter(p,'Order',1,@(x) x>=1 && isscalar(x) && mod(x,1)==0);
addParameter(p,'SeasonalityOrder',1,@(x) (isnumeric(x) && mod(x,1) == 0 && x>0) || (ischar(x) && strcmpi(x,'ByOrder')));
addParameter(p,'SeasonalityC',[]); 
addParameter(p,'SeasonalityD',[]); 
addParameter(p,'MdlProp',struct());     
addParameter(p,'IdOpt',struct());    
addParameter(p,'InitFromPreviousN',false);
addParameter(p,'FactorizeG',false,@(x) x==1 || x==0);
addParameter(p,'MinRealTol',[],@(x) isempty(x) || (isscalar(x) && isreal(x) && x>=0));
%Undoc:
addParameter(p,'Hp',[],@(x) x>0 && mod(x,1)==0 && isreal(x));
addParameter(p,'InitModel',[],@(x) isa(x,'idModels.NsfPolyModel'));

% VALIDATION
addParameter(p,'HpVal',10);                                                 
addParameter(p,'UseMae',false,@(x) x==1 || x==0);                       
% Undoc:
addParameter(p,'Sparse',false,@(x) x==0 || x==1);
addParameter(p,'ObserverInitSamples',Inf);
addParameter(p,'ObserverInitUseGls',false);
addParameter(p,'KFoldCv',1,@(x) isscalar(x) && mod(x,1)==0 && x>=1);	% Undocumented: Uses K-Fold Crossvalidation for Validation if > 1. 
parse(p,varargin{:});

%% VARIABLES
persistent M0;  % Saves Model from previous iteration.
tol = 1e-4;     % tolerance for removeing poles/zeros in A/B/E that are close to 0. 
%tol = 0;        % tolerance for removeing poles/zeros in A/B/E that are close to 0. 


%% INIT
if ~iscell(y); y = {y}; end
if isempty(u); u = cell(length(y),1); end
if ~iscell(u); u = {u}; end
Ns = length(y);
K = p.Results.KFoldCv;
n = p.Results.Order;
type = lower(p.Results.Type);
mdl_params = util.struct2namevaluepairs(p.Results.MdlProp);
opt = p.Results.IdOpt;
if ~isempty(p.Results.Hp)
   opt.Hp = p.Results.Hp; 
end
opt = util.struct2namevaluepairs(opt);


%% CHECK
assert(all(p.Results.SeasonalityC>n) && all(p.Results.SeasonalityD>n),'Seasonality needs to be > Order!');
assert(size(y{1},2)==1,'Only models with one output are supported currently!');

%% MODEL SETUP
if isfield(p.Results.MdlProp,'InputNonlinearity') && ~isempty(p.Results.MdlProp.InputNonlinearity)
    nIns = max(size(p.Results.MdlProp.InputNonlinearity));
else
    nIns = size(u{1},2);
end

switch type
    case 'arx'
        nA = n; nE = 0; nC = 0; nD = 0;
    case 'armax'
        nA = n; nE = 0; nC = n; nD = 0;
    case 'bj'
        nA = 0; nE = n; nC = n; nD = n;
    case 'oe'
        nA = 0; nE = n; nC = 0; nD = 0;
end
if ischar(p.Results.SeasonalityOrder) && strcmpi(p.Results.SeasonalityOrder,'ByOrder')
    nS = n;
else
    nS = p.Results.SeasonalityOrder;
end

ixD = false(1,max([nD+1 p.Results.SeasonalityD+1+(nS-1)]));
if strcmpi(type,'bj'); ixD(2:n+1) = true; end
for ns = 1:length(p.Results.SeasonalityD)
    ixD(p.Results.SeasonalityD(ns) + 1 + (-nS+1:nS-1)) = true;
end
ixC = false(1,max([nC+1 p.Results.SeasonalityC+1+(nS-1)]));
if any(strcmpi(type,{'armax' 'bj'})); ixC(2:n+1) = true; end
for ns = 1:length(p.Results.SeasonalityC)
    ixC(p.Results.SeasonalityC(ns) + 1 + (-nS+1:nS-1)) = true;
end

nC = length(ixC)-1;
nD = length(ixD)-1;
nB = n*ones(1,nIns);
nK = ones(1,length(nB));
M = idModels.NsfPolyModel(nA,nB,nK,nC,nD,nE,mdl_params{:});

if ~isequal(p.Results.SeasonalityD,0)
	M.D.val([false ~ixD(2:end)]) = 0; 
	M.D.free(~ixD) = false; 
end
if ~isequal(p.Results.SeasonalityC,0)
    M.C.val([false ~ixC(2:end)]) = 0; 
    M.C.free(~ixC) = false; 
end

%% INIT MODEL
if p.Results.FactorizeG
    warning off;
    M.factorize({'A' 'B' 'E'}); 
    warning on;
end

if p.Results.InitFromPreviousN && ~isempty(M0)
    try
        initParameters(M,M0,'UseRandInit',tol); 
    catch
        disp('Couldnt initialize from previous Iteration! Using specified init Method!'),
    end
end

if ~isempty(p.Results.InitModel)
    initParameters(M,p.Results.InitModel,'UseRandInit',tol);
end
    
%% IDENTIFY AND VALIDATE
for nv = 1:K
    %% SELECT DATASETS FOR VALIDATION
    ix_val = nv:K:Ns;
    if K == 1; ix_est = ix_val; else; ix_est = setdiff(1:Ns,ix_val); end
      
    %% IDENTIFICATION
    M.identify(y(ix_est),u(ix_est),opt{:});
	if p.Results.InitFromPreviousN % Save identified Model
        M0 = idModels.NsfPolyModel(nA,nB,nK,nC,nD,nE,mdl_params{:});
        if p.Results.FactorizeG
            warning off;
            M0.factorize({'A' 'B' 'E'}); 
            warning on;
        end
        initParameters(M0,M,'UseRandInit',tol);
 	end
    
    if p.Results.FactorizeG % Remove Poles in G that are close to 0 
        M.UpdateFlag = false;
        M.A = remZero(M.A,tol);
        M.B = remZero(M.B,tol);
        M.E = remZero(M.E,tol);
        M.UpdateFlag = true;
    end
    
  	if p.Results.Sparse
        M.updateSs('Sparse',true);
    end
    
    if ~isempty(p.Results.MinRealTol)
        M.updateSs('minreal',p.Results.MinRealTol);
    end
    
    %% VALIDATE
    e = M.calcResiduals(y(ix_val),u(ix_val),p.Results.HpVal,p.Results.ObserverInitSamples,~p.Results.ObserverInitUseGls); % InitMethod of observer = Inf; no GLS if true 
    E = cell2mat(e);
    
    E_k = E(:,1,end); E_k = E_k(~isnan(E_k)); % k-Step Errors
    E_1 = E(:,1,1);  E_1 = E_1(~isnan(E_1)); % 1-Step Error
    E_m = E(:); E_m = E_m(~isnan(E_m)); % k-Multistep Error
    W(nv,:) = [length(E_m) length(E_k) length(E_1)]; % Weights
    
    % Calc MMAE/MMSE of [Multistep k-Step 1-Step] predictions 
    if p.Results.UseMae
        Cost(nv,:) = [mean(abs(E_m)) mean(abs(E_k)) mean(abs(E_1))]; % Cost
    else
        Cost(nv,:) = sqrt([mean(E_m.^2) mean(E_k.^2) mean(E_1.^2)]); % Cost
    end
end
Cost = mean(Cost.*(W./sum(W,1)),1);
end

function P = remZero(P,tol)
    for r = 1:size(P,1)
        for c = 1:size(P,2)
            ixrem = find(abs(P(r,c).val(2:end))<tol)+1;  
            P(r,c).val(ixrem) = [];
            P(r,c).free(ixrem) = [];
            P(r,c).std(ixrem) = [];
        end
    end
end
    
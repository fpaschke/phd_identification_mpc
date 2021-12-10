function MDL = merge(varargin)
%MERGE Merges multiple Polynomial Models into one model.
%
% MDL = merge(mdl1,mdl2,...,mdlN)
% mdl1 ... mdlN [idModels.PolynomialModel]: Models to be merged.
% 
% EXAMPLE:
% N = 100; rng(0); U = 2*rand(N,3); Y = 2 + 1*U(:,1) + 2*U(:,2) + 3*U(:,3) + 5*randn(N,1); 
% M1 = idModels.PolynomialModel([2 2 1],'OutputName','PolyModel_1'); M1.identify(Y,U,'Nconstraints',100,'FirstDerivative',0,'SecondDerivative',0);
% M2 = idModels.PolynomialModel([2 2 1],'OutputName','PolyModel_2'); M2.identify(Y,U,'Nconstraints',100,'FirstDerivative',1,'SecondDerivative',0);
% M3 = idModels.PolynomialModel([2 2 1],'OutputName','PolyModel_3'); M3.identify(Y,U,'Nconstraints',100,'FirstDerivative',1,'SecondDerivative',1);
% M = merge(M1,M2,M3);

Nm = length(varargin);
Deg = []; InNames = {}; OutNames = {}; Incpt = []; Nv = []; InMax = []; InMin = []; OutMax = []; OutMin = [];
f = fieldnames(varargin{1}.Coefficients)'; f{2,1} = {}; Coeff = struct(f{:});
for nm = 1:Nm
    assert(isa(varargin{nm},'idModels.PolynomialModel'),'At least one model hasnt been of type "idModels.PolynomialModel"');
    assert(varargin{nm}.OutputDimension == 1,'Currently only one output models can be merged!');
    % Output related 
    OutMax = [OutMax; varargin{nm}.OutputMax];
    OutMin = [OutMin; varargin{nm}.OutputMin];
    OutNames = [OutNames; varargin{nm}.OutputName];  
    Incpt = [Incpt; varargin{nm}.OutputOffset];
    Nv = [Nv; varargin{nm}.NoiseVariance];
    % Input related
    thisInNames = varargin{nm}.InputName;
    for ni = 1:length(thisInNames)
        %check weather input already exists
        ix = strcmp(InNames,thisInNames{ni});
        if any(ix) %If Input existes 
           	InMax(nm,ix) = varargin{nm}.InputMax(ni);
            InMin(nm,ix) = varargin{nm}.InputMin(ni);
            Coeff(nm,ix) = varargin{nm}.Coefficients(1,ni); % Only Models with one Outputs supported up to now
            Deg(nm,ix) = varargin{nm}.Degree(ni);
        else % Add input
            InNames = [InNames; varargin{nm}.InputName(ni)];  
            InMax(nm,size(InMax,2)+1) = varargin{nm}.InputMax(ni);
            InMin(nm,size(InMin,2)+1) = varargin{nm}.InputMin(ni);
            Coeff(nm,size(Coeff,2)+1) = varargin{nm}.Coefficients(1,ni); % Only Models with one Outputs supported up to now
            Deg(nm,size(Deg,2)+1) = varargin{nm}.Degree(ni);
        end
    end
end
MDL = idModels.PolynomialModel(Deg,'InputName',InNames,'InputMin',min(InMin),'InputMax',min(InMax),...
                                'OutputName',OutNames,'OutputMin',OutMin,'OutputMax',OutMax,'OutputOffset',Incpt,...
                                'Coefficients',Coeff,'NoiseVariance',Nv);
end


function x0 = calcIc(obj,y,u,varargin)
%CALCIC Estimates Initial Conditions of StateSpace Representation of the LTI part of the 
%NsfPolyModel using IO Data u and y. Function preevaluates the InputNonlinearity and the uses 
%(generlized) least squares and IO-Data u and y to estimate x0.   
% 
% x0 = calcIc(obj,y,u,varargin)
% obj [NsfPolyModel]:               The Model for which Ic will be estimated.
% y [Nset cell of N x ny double]:   The Measured Outputs of the Model. Each column represents one output.
% u [Nset cell of N x nu double]:   The Measured Inputs. Each column represents one input.
% x0 [n x Nset double]:     Initial conditions at sampling instant 1. (Computation of x(0) would require u(0))
% varargin: name Value Pairs
%   'Nmax' [opt. integer]:	Number of samples from which the Ic's will be estimated. 
%                           Inf uses all available samples. Def.: 10*n
%   'UseGls' [opt. log.]:	If true then the initial conditions will be calculated by GLS. 
%                           GLS is theoratically more accurate then oridinary ls. 
%                           However the difference in most cases is not significant. Def. false
%   'FixedStates' [n x 1 double]: Can be used to fix specific (known) states.
%                           Set indices of states that should be estimated to NaN. 
%                           Example: If system has 4 states and state 1
%                           and 4 are known to be 5 and 7 then use: [5 NaN NaN 7]

%% INIT
nx = length(obj.Ss.A);
p = inputParser();
addParameter(p,'Nmax',10*nx);
addParameter(p,'UseGls',false,@(x) x == 1 || x == 0)
addParameter(p,'FixedStates',NaN(nx,1),@(x) length(x) == nx);
parse(p,varargin{:});
if strcmpi(p.Results.Nmax,'auto') || isempty(p.Results.Nmax)
    Nmax = 10*nx;
else
    Nmax = p.Results.Nmax;
end


assert(~isempty(obj.Ss),'System isnt fully parameterized which means that any of the polynomials A,B,... or E contains NaNs!');
if obj.InputDimension == 0 && (nargin <= 2 || isempty(u))
   u = zeros(length(y),0);
end
[y,u] = getRawData(obj,y,u);
u = obj.evalInputNonlinearity(u,y); % Eval inputnonlin.
y = cellfun(@(yi) bsxfun(@minus,yi,obj.OutputOffset'),y,'UniformOutput',0); %Remove OutputOffset for correct estimation

%% ESTIMATE
for ns = 1:length(y)
    if size(y{ns},1)>Nmax
        y{ns} = y{ns}(1:Nmax,:);
        if size(obj.B,2) > 0
            u{ns} = u{ns}(1:Nmax,:);
        end
    end
end
x0 = idModels.alg.lti.estIc(y,u,obj.Ss,p.Results.UseGls,obj,p.Results.FixedStates);
end
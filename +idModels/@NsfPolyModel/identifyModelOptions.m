function opt = identifyModelOptions(obj,varargin)
%IDENTIFYMODELOPTIONS Generates options struct for identify function of NsfPolyModel.
%
% opt = identifyModelOptions(obj,varargin)
% varargin: Name Value Pairs
%   'InitMethod' [char]: Determines initialization method: [Default: 'arx']
%                   'arx': uses least squares (all other coeff will be set to 0)
%                   'iv': uses instrumental variable method and PEM for Noise model
%                   'iv0': uses instrumental variable method for G and H=1
%                   'plr': uses PLR algorithm 
%                   'auto': uses iv (except if arx model is identified)
%   'InitIter' [pos int]: number of iterations to compute initial guess for parametervector if 'iv', 'plr' or 'auto' is used [Default: 20]
%   'IntegrateNoise' [ny x 1 logical]:  This Option will set an Prefilter W = 1 - q^-1 (Differnetiator) which is similar
%                                       to a Integrator in a the NoiseModel H. Thus one zero of D will be forced to 1. [Default: ny x 1 false]                                  
%   'EstimateOutputOffset' [ny x 1 logical]:    Introduces one extra parameter for the estimation of 
%                                               constant Offset of Model output y. [Default: ny x 1 false]                                  
%   'Hp' [ny x 1 pos. int]: Prediction horizon during PEM in #of samples. [Def. 1]
%   'MultiStep' [ny x 1 logical]: Uses Hp-MultiStep Criterion if true. [Def. false]
%   'Ic' [string]: Option specifying handling of initial values of residual filter. [Default: backcast]
%                   'zero': sets inital values to zero.
%                   'ls': uses least squares estimate initial values
%                   'backcast': uses backcasting to obtain initial values
%   'LsInitSamples': [pos. int. / Inf / 'auto']: Number of samples to comp. init cond. of residual filter if Ic='ls' is used. [Def. 'auto']
%                                               Inf: Uses all avail. Samples
%                                               'auto': Uses 10 times number of ic to be estimated.
%                                               pos. int.: Uses fixed number samples.
%   'LsInitWrn' [logical]: Suppress warnings which may occur during estimation of ICs with LS Method [Default: true]
%   'StabilizationMethod' [string]: Type of stabilization of optimal predictor (C,E,F Polynomials). [Def. 'none']
%                           'hard': if unstable C, E or F occurs it will be stabilized by deviding unstable modes by their
%                                   magnitude*StabilizationFactor
%                           'none': do nothing to stabilize C,E,F
%                           'constraint':   add constraints to parametervector to stabilize noisemodel. 
%   'StabilizationFactor' [scalar double >1]: Only if StabilizationMethod='hard'. See StabilizationMethod for explanation. (Def. 1.1)
%   'Solver' [string]:  Either 'matlab_***' or 'opti_***'. [Def.: 'matlab_lsqnonlin']
%                       Solver to use for solution of the problem. 
%                       matlab_***  uses solvers from optimization toolbox. Eg matlab_lsqnonlin, matlab_fmincon
%                       opti_***:   uses solvers supplied by opti toolbox. Eg opti_levmar, opti_ipopt ...
%                                   see: https://www.inverseproblem.co.nz/OPTI/%  
%   'UseGradient' ['on'/'off']:    	Use analytical gradient. [Default: 'on']
%   'CheckGradient' ['on'/'off']:   Check analytical gradient. [Default: 'off']
%                                   NOTE:   It is common that there are small errors for parameters of C/E/F at the start at each dataset. 
%                                   These errors should decay/vanish after some samples of the dataset.
%   'MaxIter' [pos int]:            Max number of iterations of optimization [Default: 500]
%   'MaxFunEval' [pos int]:         Max number of Cost function evaluations. [Default: 10e3]
%   'Display' [string]:             Show Solver output if 'iter'. See optimset() resp. optiset() for further info. [Default: 'iter']
%   'FunTolAbs' [pos. double]:      Abortion criterion of optimization. Abort if |F_i - F_i-1| <= FunTolAbs [Default: 1e-8]
%   'TolX'  [pos, double]:          Abortion criterion of optimization.  Stops if norm of change in parametervector is less
%                                   then TolX [Default: 1e-8]
%   'ForcePolesG' [char]:           Forces poles of deterministic subsystem G(q) to be
%                                   's': stable 
%                                	'p': positive real 
%                                   'sp' or 'ps': positive real and stable
%                                   '': No constraints on poles/zeros [Default]
%                                   by introducing appropriate constraints.
%   'ForceSymmetricA' [logical]: 	Uses equality constraints to ensure symmetric A Matrix                              
%   NOTE: ALL UNMATCHED PARAMETERS WILL BE PASSED TO optimset() (matlab solver) OR TO optiset() (opti solver)
%
%   UNDOCUMNETED (LEAVE UNCHANGED):
%   'Feedback' ['raw' or 'sim']:    In case of models with outputfeedback the inputnonlinearity will be evaluated as [Def. 'raw']:
%                                   'raw' f(u[t],y[t]) where y is the measured data of the output 
%                                   'sim' f(u[t],y_hat[t|t-Hp]) where y_hat will be computed by calcPredictions

p = inputParser(); p.KeepUnmatched = true;
addParameter(p,'InitMethod','arx',@(x) any(strcmpi(x,{'plr' 'iv' 'iv0' 'arx' 'auto'}))); 
addParameter(p,'InitIter',20,@(x) mod(x,1)==0 && x>=1);
addParameter(p,'IntegrateNoise',false(obj.OutputDimension,1),@(x) all(x==1 | x == 0) && length(x)==obj.OutputDimension || isscalar(x));
addParameter(p,'Hp',ones(obj.OutputDimension,1),@(x) all(mod(x,1)==0 & x>=1) && (length(x) == obj.OutputDimension) || isscalar(x));
addParameter(p,'MultiStep',false,@(x) x==1 || x == 0);
addParameter(p,'Ic','backcast',@(x) any(strcmpi(x,{'zero' 'ls' 'backcast'})));
addParameter(p,'LsInitWrn',false,@(x) x==1 || x==0);
addParameter(p,'LsInitSamples','auto',@(x) strcmpi(x,'auto') || (isnumeric(x) && (isinf(x) || mod(x,1)==0)));
addParameter(p,'StabilizationMethod','none',@(x) any(strcmpi(x,{'none' 'constraint' 'hard'}))) % HARD: Hat nen BUG, bei teilweise festem C
addParameter(p,'StabilizationFactor',1.1,@(x) x>1)
addParameter(p,'Solver','matlab_lsqnonlin',@ischar);
addParameter(p,'CheckGradient','off',@(x) any(strcmpi(x,{'on' 'off'})));
addParameter(p,'UseGradient','on',@(x) any(strcmpi(x,{'on' 'off'})));
addParameter(p,'Display','iter',@ischar);
addParameter(p,'MaxIter',5e2,@(x) (mod(x,1)==0 && x>1) || isinf(x));
addParameter(p,'MaxFunEval',Inf,@(x) (mod(x,1)==0 && x>1) || isinf(x))
addParameter(p,'FunTolAbs',1e-8,@(x) isscalar(x) && isnumeric(x) && x>0);
addParameter(p,'TolX',1e-8,@(x) isscalar(x) && isnumeric(x) && x>0);
addParameter(p,'ForcePolesG','',@(x) any(strcmpi(x,{'' 'p' 's' 'ps' 'sp'})));
addParameter(p,'ForceZerosG','',@(x) any(strcmpi(x,{'' 'p' 's' 'ps' 'sp'}))); 
addParameter(p,'ForceSymmetricA',false,@(x) x==1 || x==0); 

%
addParameter(p,'PlotZeros',{},@(x) iscell(x) || any(strcmpi(x,{'A' 'C' 'D' 'E'})))
addParameter(p,'Feedback','raw',@(x) any(strcmpi(x,{'raw' 'sim'})));
parse(p,varargin{:});
opt = p.Results;
opt.Unmatched = p.Unmatched;

if isscalar(opt.Hp); opt.Hp = opt.Hp*ones(obj.OutputDimension,1); opt.Hp = opt.Hp(:); end 
if isscalar(opt.IntegrateNoise); opt.IntegrateNoise = logical(opt.IntegrateNoise*ones(obj.OutputDimension,1)); end 
if ischar(opt.PlotZeros); opt.PlotZeros = {opt.PlotZeros}; end
end
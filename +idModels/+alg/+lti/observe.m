function [y_hat, x_hat, e] = observe(S,Y,U,X0,y0)
%OBSERVE Runs Observer for LTI system in innovations form.
% x[k+1] = Ax_hat[k] + Bu[k] + K(y[k]-y_hat[k])
% y_hat[k]   = Cx_hat[k] + Du[k]
%where the gain of the observer K is given by the noise part of the system S. 
% 
% [y_hat, x_hat, e] = observe(S,Y,U,X0,y0)
% S [ss]:   Matlab Statespace object. Requires S to be in innovations form, which means that S.InputGroup.Noise needs to be 
%           available. See matlab docoumentation of ss() for further info.
% Y [Nset cell of N x ny double]: The Measured Outputs of the Model. Each column represents one output. Y = [y(1); ...; y(N)]
% U [Nset cell of N x nu double]: The Measured Inputs. Each column represents one input. U = [u(1); ...; u(N)]
% X0 [nx x Nset double]: Matrix of initial conditions at sampling instant ("index") 1.
% y0 [opt. ny x 1 double]: Outputoffset of ss Model. Def.: zeros(ny,1). 
% y_hat [Nset cell of N x ny double]: The Observed Output of the Model. y_hat = [y_hat(1); ...; y_hat(N)]
% x_hat [Nset cell of N x nx double]: The Observed state of the Model. y_hat = [x_hat(1)=X0; ...; x_hat(N)]
% e [Nset cell of N x nx double]: The Obervation Error e: e = y - y_hat.

%% PREPARE DATA AND INIT
ny = size(S.C,1);
if nargin<6 || isempty(DoChecks)
    DoChecks = true;
end
if nargin == 4
    y0 = zeros(size(S.C,1),1);
else
    assert(isvector(y0),'y0 needs to be a double vector!');
    y0 = y0(:);
end
if ~iscell(Y)
   Y = {Y}; 
end
Nsets = length(Y);
if ~iscell(U)
    if isempty(U)
        U = cell(Nsets,1);
    else
        U = {U}; 
    end
end

if isempty(U)
   U = cell(Nsets,1); 
end

if isfield(S.InputGroup,'Measured')
    ix_u = S.InputGroup.Measured;
else
    ix_u = [];
end

%% RUN OBSERVATION
KG = idModels.alg.lti.kalmanGain(S);
Ao = S.A-KG*S.C;
Bo = [S.B(:,ix_u)-KG*S.D(:,ix_u) KG];
Co = S.C;
Do = [S.D(:,ix_u) zeros(size(S.D,1),ny)];

for ns = 1:Nsets   
	[y_hat{ns,1}, x_hat{ns,1}] = idModels.alg.lti.lsim_ss(Ao,Bo,Co,Do,[U{ns} bsxfun(@minus,Y{ns},y0')],sparse(X0(:,ns))); %yo(1),... ; xo(1),...
	y_hat{ns,1} = bsxfun(@plus,y_hat{ns},y0');
    if nargout == 3
        e{ns,1} = Y{ns} - y_hat{ns};  
    end
end
end
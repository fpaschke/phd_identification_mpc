function [x0,PSI,THETA_e] = estIc(Y,U,sys,UseGls,mdl_poly,states)
%ESTIC Estimate x(1) from measured i/o data y(1),y(2), ... and u(1), u(2), ...; 
% x(0) cant be estimated since it would require u(0)...
%
%x0 = estIc(y,u,sys,init_oe,mdl_poly)
% Estimate of Initial Conditions assuming augmented stochastic statespace system:
% x(k+1) = A*x(k) + B*u(k) + K*e(k)
% y(k) = C*x(k) + D*u(k) + e(k)
% If the model is purely deterministic (There is no field sys.InputGroup.Noise) then standard LS is used to estimate initial
% conditions (OE-Model is assumed)
%
% y [Nxny double or cell with Nxny double]: Measurements of output each cell corresponds to one experiment. If [] then 0
%                                           will be assumed
% u [Nxnu double or cell with Nxnu double]: Measurements of input each cell corresponds to one experiment
% sys [ss]: Matlab Statespace Object.
% UseGls [opt. logical]:   If true then the initial conditions will be calculated by GLS. 
%                        	GLS is theoratically more accurate then oridinary ls. 
%                           However the difference in most cases is not significant. Def. false
% mdl_poly [opt. NsfPolyModel]:	Function uses transfer function description to simulate forced part of the system, 
%                            	which are computed from the the statespace object sys. If the model is large this can
%                            	be timeconsuming. In these cases the polynomial description of the system "mdl_poly"
%                             	can be supplied, from which the computation of the transfer functions is much faster.

%% INIT
%put to predefined scheme
if ~iscell(Y) && ~isempty(Y)
    Y = {Y}; 
end
if ~iscell(U) && ~isempty(U)
    U = {U}; 
end
if isempty(Y)
    Nexp = length(U);
    N = cellfun(@(u) size(u,1),U);
else
    Nexp = length(Y);
    N = cellfun(@(y) size(y,1),Y);
end
Nmax = max(N);
if nargin < 4 || isempty(UseGls); UseGls = false; end 

if isprop(sys,'InputGroup') && isfield(sys.InputGroup,'Noise')
    ix_e = sys.InputGroup.Noise;
    nu = size(sys.B,2) - length(ix_e);
    SigmaE = sys.D(:,ix_e);
    CovE = SigmaE*SigmaE'; 
else
    ix_e = [];
    nu = size(sys.B,2);
    UseGls = false;
end
ix_u = setdiff(1:size(sys.B,2),ix_e);
A = sys.A; B = sys.B(:,ix_u); C = sys.C; D = sys.D(:,ix_u); K = sys.B(:,ix_e);
nx = size(A,1);
ny = size(sys.C,1);    
if nargin <= 5 || isempty(states)
    states = NaN(nx,1);
end
ix_est = isnan(states);

%% CALC TF SYSTEM
PSI = idModels.alg.lti.calcSsPredictors(A,[],C,[],Nmax-1);
if nargin <= 4 % we need to comp tfs
    if nu>0
        G = tf(ss(A,B,C,D,1)); G.Variable = 'z^-1';
        G_num = G.num; G_den = G.den;
    end
    if UseGls
        H = tf(ss(A,K,C,SigmaE,1)); H.Variable = 'z^-1';
        H_num = H.num; H_den = H.den;
    end
else 
    if UseGls
        [G_num,G_den,H_num,H_den] = mdl_poly.calcGH('AbsorbNoiseVariance',false);
    else
        [G_num,G_den] = mdl_poly.calcGH('AbsorbNoiseVariance',false);
    end
end

if UseGls
    assert(isdiag(CovE));
    hd = NaN(ny,Nmax);
    for no = 1:ny
     	hd(no,:) = filter(H_num{no,no},H_den{no,no},[1; zeros(Nmax-1,1)]);
    end
    THETA_e = tril(toeplitz(hd(:)));
    R = THETA_e*THETA_e'; %covariance of Noise
end

%if rank(PSI)<min(size(PSI)); warning('PSI-Matrix is rank deficient!'); end 
%wS = warning(); warning off; % to avoid multiple warnings if PSI is rank deficient

%% EST
x0 = NaN(nx,Nexp);
for ne = 1:Nexp
    % yf = THETA_u(1:N(ne)*ny,1:N(ne)*nu)*u; %slow
    x0(~ix_est,ne) = states(~ix_est); 
    yf = zeros(ny,N(ne));
	if nu > 0  
        for no = 1:ny
            yf(no,:) = sum(cell2mat(arrayfun(@(ni) filter(G_num{no,ni},G_den{no,ni},U{ne}(:,ni)),1:nu,'UniformOutput',0)),2)';
        end
 	end
    yf = yf(:);
    if ~isempty(Y)
        y = Y{ne}'; y = y(:);
        ed = y - yf;
    else
        ed = -yf;
    end
    if any(~ix_est)
       ed = ed - PSI(1:N(ne)*ny,~ix_est)*states(~ix_est); 
    end
    
    if ~UseGls     
        x0(ix_est,ne) = PSI(1:N(ne)*ny,ix_est)\ed;
    else
        x0(ix_est,ne) = lscov(PSI(1:N(ne)*ny,ix_est),ed,R(1:N(ne)*ny,1:N(ne)*ny));
    end
end
%warning(wS); %restore warning state
end
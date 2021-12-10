function [PSI, THETA_U, THETA_E] = calcSsPredictors(A,B,C,D,Hp,K)
%Computes Matrices PSI and THETA_U and THETA_E for calculation of Multistepprediction 
% Y = PSI*x0 + THETA_U*U
%of a statespace Model 
% x(k+1) = A*x(k) + B*u(k) + K*e(k)
% y(k) = C*x(k) + D*u(k) + e(k)
%defined by A B C D and K where 
% U = [u(k); u(k+1); ... ; u(k+Hp)] and Y 
% Y = [y(k); ... ; y(k+Hp)] are the vectors of input resp. outputpredictions.
% PSI = [C, CA, ... CA^Hp]'; 
% THETA_U = [D 0 ...  0; CB D 0 ...0; CAB CB D 0 ... 0; ... ; CA^Hp-1B ...]'
% THETA_E = [I 0 ...  0; CK I 0 ...0; CAK CK I 0 ... 0; ... ; CA^Hp-1K ...]'
%
% [PSI, THETA_U, THETA_E] = calcSsPredictors(A,B,C,D,Hp,K)
% A,B,C,D,K [double Matrices]: The Statespace matrices of the system. 
% Hp [pos. integer]: Desired Prediction Horizon in number of steps.
% PSI [(Hp+1)*ny x nx double]: Extended observability Matrix.
% THETA_U [ny*(Hp+1) x nu*(Hp+1)]: Toplitz Matrix Containing the Markov Parameters of the deterministic part of the system.
% THETA_E [ny*(Hp+1) x ny*(Hp+1)]: Toplitz Matrix Containing the Markov Parameters of the stochastic part of the system.
%
% [PSI, THETA_U, THETA_E] = calcSsPredictors(sys,Hp)
% sys [ss]: The Statespace description of the system. 
% Hp [pos. integer]: Desired Prediction Horizon in number of steps.
% PSI [(Hp+1)*ny x nx double]: Extended observability Matrix.
% THETA_U [ny*(Hp+1) x nu*(Hp+1)]: Toplitz Matrix Containing the Markov Parameters of the deterministic part of the system.
% THETA_E [ny*(Hp+1) x ny*(Hp+1)]: Toplitz Matrix Containing the Markov Parameters of the stochastic part of the system.

if nargin == 2
    Hp = B;
    sys = A;
    if isfield(sys.InputGroup,'Noise')
        ix_e = sys.InputGroup.Noise;
    else
        ix_e = [];
    end
    ix_u = setdiff(1:size(sys.B,2),ix_e);
  	A = sys.A; B = sys.B(:,ix_u); C = sys.C; D = sys.D(:,ix_u); K = sys.B(:,ix_e);
end

ny = size(C,1);
nu = size(B,2);
n = size(A,1);

if isempty(D)
    D = zeros(ny,nu);
end

PSI = zeros(ny*(Hp+1),n);
if nargout > 1
    THETA_U = zeros(ny*(Hp+1),nu*(Hp+1));
    THETA_U(1:ny,1:nu) = D;
end
if nargout > 2
    THETA_E = zeros(ny*(Hp+1),ny*(Hp+1));
    THETA_E(1:ny,1:ny) = eye(ny);
end

Ai = A; Aim1 = sparse(eye(n));
PSI(1:ny,:) = C;

for i = 1:Hp
    rows = i*ny+1:(i+1)*ny;
    %PSI(rows,:) = C*Ai;
    PSI(rows,:) = PSI(rows-ny,:)*A;
    if nargout > 1
        THETA_U(rows,1:(i+1)*nu) = [C*Aim1*B THETA_U(rows-ny,1:i*nu)];
      	Aim1 = Ai;
        Ai = Ai*A;
    end
    if nargout > 2
      	THETA_E(rows,1:(i+1)*ny) = [C*Aim1*K THETA_E(rows-ny,1:i*ny)];
    end
end
end
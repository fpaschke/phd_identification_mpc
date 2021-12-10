clc; clear; close all;

%% PARAMS AND DATA
% Params
nu = 2;
A = {[-.95 .5; 0 -.9]}; ny = size(A{1},2);
B = {zeros(ny,nu) eye(ny,nu)}; 
C = A;

% Data
Np = 10; Ns = 10;
std_e = .2;
u = kron(rand(Np,nu),ones(Ns,1));
e = std_e*randn(Np*Ns,ny);

%% SIM AND PLOT
y_ = idModels.alg.lti.lsim_poly(A,B,[],u);
y = y_ + e;

subplot(2,1,1); plot(y(:,1),'-k'); hold on; plot(y_(:,1),'--k'); yyaxis right;
plot(y(:,2),'-b'); hold on; plot(y_(:,2),'--b');
subplot(2,1,2); plot(u);

%% IDENT
[A,B,var_e,e] = idModels.alg.ls.oe_iv(y,u,ones(ny),eye(ny,nu),ones(ny,nu),'Plot',true);
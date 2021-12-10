clc; clear; close all;

%% PARAMS AND DATA
% Params
A = {[-.9 .5; 0 -.8]}; ny = size(A{1},2);
B = {zeros(ny,2) eye(ny,2)}; nu = size(B{1},2);
C = {[.8 .0; .0 .7]};

% Data
rng(0)
Np = 20; Ns = 50;
std_e = [.1 .1];
u = kron(rand(Np,nu),ones(Ns,1));
e = std_e.*randn(Np*Ns,ny);

%% SIM AND PLOT
y_ = idModels.alg.lti.lsim_poly(A,B,C,u);
y = idModels.alg.lti.lsim_poly(A,B,C,u,e);

subplot(2,1,1); plot(y(:,1),'-k'); hold on; plot(y_(:,1),'--k'); yyaxis right;
plot(y(:,2),'-b'); hold on; plot(y_(:,2),'--b');
subplot(2,1,2); plot(u);

%% IDENT
[A,B,C,D,E,var_e] = idModels.alg.ls.plr(y,u,ones(ny),eye(ny,nu),ones(ny,nu),eye(ny),'Iter',50,'Plot',true);
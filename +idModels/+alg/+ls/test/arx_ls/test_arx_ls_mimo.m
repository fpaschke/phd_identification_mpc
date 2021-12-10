clc; clear; close all;

%% PARAMS AND DATA
% Params
nu = 2;
A = {[-.95 -.1; 0 -.9]}; ny = size(A{1},2);
B = {zeros(ny,nu) eye(ny,nu)}; 
C = {eye(ny)};

% Data
Np = 50; Ns = 20;
std_e = .1;
u = kron(rand(Np,nu),ones(Ns,1));
e = std_e*randn(Np*Ns,ny);

%% SIM AND PLOT
y_ = idModels.alg.lti.lsim_poly(A,B,[],u);
y = idModels.alg.lti.lsim_poly(A,B,[],u,e);

subplot(2,1,1); plot(y(:,1),'-k'); hold on; plot(y_(:,1),'--k'); yyaxis right;
plot(y(:,2),'-b'); hold on; plot(y_(:,2),'--b');
subplot(2,1,2); plot(u);

%% IDENT
a = idModels.util.polyStruct([1 1; 1 1],'A');
b = idModels.util.polyStruct([1 1; 1 1],'B',[1 1; 1 1]);
[A_hat,B_hat,var_e,e] = idModels.alg.ls.arx_ls(y,u,a,b);
A_hat{1,:}
A_hat{2,:}
B_hat{1,:}
B_hat{2,:}
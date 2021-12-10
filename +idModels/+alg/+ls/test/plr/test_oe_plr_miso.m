clc; clear; close all;

%% PARAMS AND DATA
% Params
nu = 2; ny = 1;
B11 = [0 1]; B12 = B11;
F1 = [1 -.9];

% Data
rng(0)
Np = 20; Ns = 50;
std_e = [.1];
u = kron(rand(Np,nu),ones(Ns,1));
e = std_e.*randn(Np*Ns,ny);

%% SIM AND PLOT
y_ = filter(B11,F1,u(:,1)) + filter(B12,F1,u(:,2));
y = y_ + e;

subplot(2,1,1); plot(y(:,1),'-k'); hold on; plot(y_(:,1),'--k'); 
subplot(2,1,2); plot(u);

%% IDENT OE
[A,B,C,D,E,var_e] = idModels.alg.ls.plr(y,u,[],ones(ny,nu),ones(ny,nu),[],[],eye(ny),'Iter',50,'Plot',true); %-> divergent 

%% IDENT ARMAX
[A,B,C,D,E,var_e] = idModels.alg.ls.plr(y,u,1,ones(ny,nu),ones(ny,nu),1,[],[],'Iter',50,'Plot',true); %-> works
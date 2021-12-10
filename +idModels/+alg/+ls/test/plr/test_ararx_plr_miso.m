clc; clear; close all;

%% PARAMS AND DATA
% Params
nu = 2; ny = 1;
A = conv([1 -.9],1);
B1 = [0 1]; B2 = [0 .5]; 
D = [1 -.6];
AD = conv(A,D);

% Data
rng(0)
Np = 1000; Ns = 1;
std_e = .1;
u = kron(rand(Np,nu),ones(Ns,1));
e = std_e*randn(Np*Ns,ny);

%% SIM AND PLOT
y_ = filter(B1,A,u(:,1)) + filter(B2,A,u(:,2));
y = y_ + filter(1,AD,e);

subplot(2,1,1); plot(y_(:,1),'-k'); hold on; plot(y(:,1),'--r');
subplot(2,1,2); plot(u);

%% IDENT
[A,B,C,D,F,var_e,e] = idModels.alg.ls.plr(y,u,1,[1 1],[1 1],0,1,0,'plot',1);
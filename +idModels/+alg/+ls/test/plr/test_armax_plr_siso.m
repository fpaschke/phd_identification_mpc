clc; clear; close all;

%% PARAMS AND DATA
% Params
A = conv([1 -.8],[1 -.7]) 
B = [0 1];
C = [1 .7];

% Data
rng(0)
Np = 10; Ns = 100;
std_e = .01;
u = kron(rand(Np,1),ones(Ns,1));
e = std_e*randn(Np*Ns,1);

%% SIM AND PLOT
y_ = filter(B,A,u(:,1));
y = y_ + filter(C,A,e);

subplot(2,1,1); plot(y_(:,1),'-k'); hold on; plot(y(:,1),'--r');
subplot(2,1,2); plot(u);

%% IDENT ARMAX
[A,B,C,D,F,var_e,e] = idModels.alg.ls.plr(y,u,2,1,1,1,'plot',1);
%% IDENT ARIMAX
[A,B,C,D,F,var_e,e] = idModels.alg.ls.plr(y+20,u,2,1,1,1,'plot',1,'Prefilter',{[1 -1] 1});
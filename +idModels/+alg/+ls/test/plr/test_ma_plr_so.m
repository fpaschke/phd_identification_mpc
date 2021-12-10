clc; clear; close all;

%% PARAMS AND DATA
% Params
C = conv([1 -.9],[1 0]); 

% Data
rng(0)
e = .1*randn(1e3,1);

%% SIM AND PLOT
y = filter(C,1,e);
plot(y(:,1),'-k');

%% IDENT
[A_,B_,C_] = idModels.alg.ls.plr(y,[],0,[],[],2,'Iter',50,'Plot',true);
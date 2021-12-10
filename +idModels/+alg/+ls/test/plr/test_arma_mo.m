clc; clear; close all;

%% PARAMS AND DATA
% Params
nu = 2;
A = {[-.9 .5; 0 -.8]}; ny = size(A{1},2);
C = {[.8 .0; .0 .7]};

% Data
rng(0)
Np = 20; Ns = 50;
std_e = [.1 .1];
e = std_e.*randn(Np*Ns,ny);

%% SIM AND PLOT
y = idModels.alg.lti.lsim_poly(A,[],C,[],e);
plot(y(:,1),'-k');  yyaxis right; plot(y(:,2),'-b');

%% IDENT
[A,B,C,var_e] = idModels.alg.ls.plr(y,[],ones(ny),[],[],eye(ny),'Iter',50,'Plot',true);
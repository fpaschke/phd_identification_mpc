clc; clear; close all;

%% PARAMS
F = [1 -.95]; 
B = [0 .1];
rng(0); u = kron(rand(10,1),ones(10,1));
yd = filter(B,F,u,1);
y = yd + .1*randn(length(u),1);
plot([yd y(:,1)])

%% ESTIMATE OE
[Ai,Bi] = idModels.alg.ls.oe_iv(y,u,1,1,1,'EstIc',true,'Plot',true);

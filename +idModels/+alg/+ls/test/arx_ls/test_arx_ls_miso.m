clc; clear; close all;

%% PARAMS AND DATA
% Params
nu = 2; ny = 1;
A = conv([1 -.8],[1 -.7])
B = {[0 1] [1 0]}; 

% Data
Np = 100; Ns = 100;
std_e = .1;
u = kron(rand(Np,nu),ones(Ns,1));
e = std_e*randn(Np*Ns,ny);

%% SIM AND PLOT
y_ = filter(B{1},A,u(:,1)) + filter(B{2},A,u(:,2));
y = y_ + filter(1,A,e);
 
subplot(2,1,1); plot(y(:,1),'-k'); hold on; plot(y_(:,1),'--r');
subplot(2,1,2); plot(u);

%% IDENT ARX
a = idModels.util.polyStruct(2,'A');
b = idModels.util.polyStruct([1 1],'B',[1 0]); b(2).val(2) = 0; b(2).free(2) = false;
[A_arx,B_arx,var_e,e] = idModels.alg.ls.arx_ls(y,u,a,b);
% [A_arx,B_arx,var_e,e] = idModels.alg.ls.arx_ls(y,u,2,[1 0],[1 0]); %equivalent
A_arx{1}
B_arx{1}
B_arx{2}

%% IDENT ARIX
[A_arix,B_arix,var_e,e] = idModels.alg.ls.arx_ls(y+10,u,a,b,[],'PreFilter',{[1 -1] 1});
A_arix{1}
B_arix{1}
B_arix{2}
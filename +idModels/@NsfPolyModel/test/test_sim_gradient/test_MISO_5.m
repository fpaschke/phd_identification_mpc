clc;
import idModels.*

%% MISO Nonlinear Model
M = load('TestModel.mat'); M = M.Mod;
rng(0) %seed random number generator

% generate some data for testing
m_Tsup = kron([rand(10,1) 5*rand(10,1)+18],ones(20,1)); 
N = length(m_Tsup);
u = [m_Tsup [5 0 50 50].*ones(N,1)];

% setup model
u_grads = [1 2]; % Compute Gradient wrt u1 and u2
[y,x,dy_du] = M.simulate(u,'CalcInputGradients',u_grads); 

%% Check with finite differences
figure('color','w'); 
j = 1; % 
delta = 1e-6;
i=1;
for k = u_grads
    up = u; um = u;
    up(j,k) = up(j,k) + delta; 
    um(j,k) = um(j,k) - delta; 
    grad = (M.simulate(up) - M.simulate(um))/(2*delta); % numerical gradient
    Dif = grad -squeeze(dy_du(:,1,i,j)); % shoul be very small (>1e-6)grad - 
    subplot(size(u_grads,2),1,i); plot(Dif); ylabel(['\Delta dy/du_{' num2str(k) '}(' num2str(j) ')'])
    yyaxis right; plot(squeeze(dy_du(:,1,i,j))); ylabel(['dy/du_{' num2str(k) '}(' num2str(j) ')'])
    i=i+1;
end
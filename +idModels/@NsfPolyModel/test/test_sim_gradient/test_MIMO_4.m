clc;
import idModels.*

%% MISO Nonlinear Model
load('TestModel3.mat'); 

% generate some data for testing
rng(0) %seed random number generator
u = kron([zeros(10,4) 10*rand(10,1)+20 4*rand(10,4)],ones(10,1));
u(:,6:9) = [ones(1,4); diff(u(:,6:9),1)];

% setup model
u_grads = 5:9; % Compute Gradient wrt u5 ... u9
[y,x,dy_du] = BuilMod.simulate(u,'CalcInputGradients',u_grads); 

%% Check with finite differences
figure('color','w'); 
j = 10; 
delta = 1e-3;
m=1;
ny = 4;
for l=1:ny
    i=1;
    for k = u_grads
        up = u; um = u;
        up(j,k) = up(j,k) + delta; 
        um(j,k) = um(j,k) - delta;
        grad = (BuilMod.simulate(up) - BuilMod.simulate(um))/(2*delta); % numerical gradient
        Dif = grad - squeeze(dy_du(:,:,i,j)); % shoul be very small (>1e-6)grad - 
        subplot(ny,length(u_grads),m); 
        plot(Dif(:,l)); hold on; ylabel(['\Delta dy_{' num2str(l) '}/du_{' num2str(k) '}(' num2str(j) ')']);
        yyaxis right; plot(squeeze(dy_du(:,l,i,j))); ylabel(['dy_{' num2str(l) '}/du_{' num2str(k) '}(' num2str(j) ')']);
        m = m + 1;
        i = i + 1;  
    end   
end
import idModels.*

close all; 
clc; clear;

% % % Nichtlinearitšt f(u,y) = u(1)*(p1-y) definieren
F1 = @(u,p,y) deal(u.*(p(1) - y),... %f
    u,... %df/dp
    (p(1) - y),... %df/du
    -u); %df/dy
%f2 = u2^2*y2+u2^3*p1
F2 = @(u,p,y) deal( u.^2.*y + u.^3.*p(1),...
                    u.^3,...
                    2*y.*u+3*p(1).*u.^2,...
                    u.^2);

I(1).fun = F1; I(1).parameters = 22; I(1).input_idx = [1]; I(1).output_idx = [1];
I(2).fun = F2; I(2).parameters = [15]; I(2).input_idx = [2]; I(2).output_idx = [2];

% Modell instanzieren
M = idModels.NsfPolyModel(eye(2),eye(2),eye(2),'InputNonlinearity',I,'OutputOffset',[20 10]');
M.A(1,1).val = [1 -.9]; M.B(1,1).val = [0 .1];  
M.A(2,2).val = [1 -.5]; M.B(2,2).val = [0 .5]; 
M.NoiseVariance = .1*eye(2);
u_grads=[1 2];

% Modell Simulieren
rng(0)
u = kron(rand(10,2),ones(20,1));
[y,x,dy_du] = M.simulate(u,'CalcInputGradients',u_grads);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

figure('Color','w');    
subplot(2,1,1); plot(y(:,1)); grid on;
subplot(2,1,2); plot(y(:,2)); grid on;  

%% Check with finite differences
figure('color','w'); 
j = 50; % 
delta = 1e-8;
m=1;
for l=1:size(y,2)
    i=1;
    for k = u_grads
        up = u; um = u;
        up(j,k) = up(j,k) + delta; 
        um(j,k) = um(j,k) - delta;
        grad = (M.simulate(up) - M.simulate(um))/(2*delta); % numerical gradient
        Dif = grad - squeeze(dy_du(:,:,i,j)); % shoul be very small (>1e-6)grad - 
        subplot(length(u_grads),size(y,2),m); 
        plot(Dif(:,l)); hold on; ylabel(['\Delta dy_{' num2str(l) '}/du_{' num2str(k) '}(' num2str(j) ')']);
        yyaxis right; plot(squeeze(dy_du(:,l,i,j))); ylabel(['dy_{' num2str(l) '}/du_{' num2str(k) '}(' num2str(j) ')']);
        m = m + 1;
        i = i + 1;  
    end   
end
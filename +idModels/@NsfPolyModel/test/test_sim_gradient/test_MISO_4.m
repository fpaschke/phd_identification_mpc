close all; clear; clc;
import idModels.*

%% MISO Nonlinear Model
%rng(0) %seed random number generator

% generate some data for testing
A = conv([1 -.8+.5i],[1 -.8-.5i]); na = length(A)-1; 
u1 = kron(randn(1,10),ones(1,5))'; 
u2 = kron(randn(1,10),ones(1,5))';
u3 = kron(randn(1,10),ones(1,5))';
u = [u1 u2 u3];

% Coefficients of 
A = conv([1 -.8],[1 -.7]);
B = {[0 0.1] [0 .02] [0 0.05]};
C = [1 .7];

% setup inpput Nonlinearity (Vector field f(u(k),y(k))

% f1 = u1^2*u2*y
f(1).parameters = [];
f(1).input_idx = [1 2];
f(1).output_idx = [1];
f(1).fun = @(u,p,y) deal( u(:,1).^2.*u(:,2).*y,...
                        [],...
                        [2*u(:,1).*u(:,2).*y u(:,1).^2.*y],...
                        u(:,1).^2.*u(:,2));


% f2 = p1*u2^2+p2*u3^3
f(2).parameters = [2 -1];
f(2).input_idx = [2 3];
f(2).output_idx = [];
f(2).fun = @(u,p) deal(   p(1).*u(:,1).^2 + p(2).*u(:,2).^3,...
                          [],...
                          [2*p(1).*u(:,1) 3*p(2).*u(:,2).^2 ],...  %%%%%%%%%%%%%%
                          []); 

% f3 = u1*u2 + u3*y
f(3).parameters = [];
f(3).input_idx = [1 2 3];
f(3).output_idx = 1;
f(3).fun = @(u,p,y) deal( u(:,1).*u(:,2) + u(:,3)*y,...
                        [],...
                        [u(:,2) u(:,1) y],...   %%%%%%%%%%%%%%%
                        u(:,3)); 


% setup model
M1 = idModels.NsfPolyModel(2,[1 1 1],[1 1 1],1,'InputNonlinearity',f);
M1.A.val = A; M1.B(1,1).val = B{1}; M1.B(1,2).val = B{2}; M1.B(1,3).val = B{3}; M1.C.val = C; M1.NoiseVariance = 0; % Set coefficients of polynomials
u_grads=3;
[y,x,dy_du] = M1.simulate(u,'CalcInputGradients',u_grads); % Compute Gradient wrt u1 u2 and u3 

% Plot model response
subplot(2,1,1); plot(y); ylabel('y'); xlabel('Time k'); grid on;
subplot(2,1,2); plot(u); ylabel('u'); xlabel('Time k'); legend({'u1' 'u2' 'u3'}); grid on;

%% Check with finite differences
figure('color','w'); 
j = 10; 
delta = 1e-8;
i=1;
for k = u_grads
    up = u; um = u;
    up(j,k) = up(j,k) + delta; 
    um(j,k) = um(j,k) - delta; 
    grad = (M1.simulate(up) - M1.simulate(um))/(2*delta); % numerical gradient
    Dif = grad -squeeze(dy_du(:,1,i,j)); % shoul be very small (>1e-6)grad - 
    subplot(length(u_grads),1,i); plot(Dif); ylabel(['\Delta dy/du_{' num2str(k) '}(' num2str(j) ')'])
    yyaxis right; plot(squeeze(dy_du(:,1,i,j))); ylabel(['dy/du_{' num2str(k) '}(' num2str(j) ')'])
    i=i+1;
end
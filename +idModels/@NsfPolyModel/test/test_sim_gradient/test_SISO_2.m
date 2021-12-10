close all; clear; clc;
import idModels.*
%import time.*
% Polynome A und B
A = conv([1 -.8],1); B = [0 .3]; y0 = 10;

na = length(A)-1; nk = find(B~=0,1,'first')-1; nb = length(B)-nk; 

% % % linear f(u,y) = u(1) definieren
F = @(u,p,y) deal(u,... %f
    [],... %df/dp
    ones(size(u,1),1),... %df/du
    zeros(size(u,1),1)); %df/dy
I(1).fun = F; 
I(1).parameters = [];
I(1).free = []; 
I(1).input_idx = 1; 
I(1).output_idx = [];
% Modell simulieren
u1 = kron(rand(10,1),ones(15,1));

x(1,1) = 0; y = y0;
for k = 2:length(u1)
     [fu,~,~,~] = F(u1(k-1:-1:k-nb),[],y(k-1:-1:k-nb));
     x(k,1) = -A(2:end)*x(k-1:-1:k-na) + B(2:end)*fu;
     y(k,1) = y0 + x(k,1);
end

% Modell identifizieren
M = idModels.NsfPolyModel(1,1,1,'InputNonlinearity',I);
M.identify(y,u1,'EstimateOutputOffset',1,'FunTolAbs',1e-20); 

% Modell Simulieren
[ys,x,dy_du] = M.simulate(u1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
plot(y); hold on; plot(ys,'--r')

%% Check with finite differences
figure('color','w'); 
j = 20; %Zeitschritt f�r u f�r welchen der gradient berechnet werden soll!
delta = 1e-8;
for k = 1:size(u1,2)
    up = u1; um = u1;
    up(j,k) = up(j,k) + delta; 
    um(j,k) = um(j,k) - delta; 
    grad = (M.simulate(up) - M.simulate(um))/(2*delta); % numerical gradient
    Dif = grad - squeeze(dy_du(:,1,k,j)); % should be very small (>1e-6) 
    subplot(size(u1,2),1,k); 
    plot(Dif);
    ylabel(['dy/du(' num2str(j) ')'])
end
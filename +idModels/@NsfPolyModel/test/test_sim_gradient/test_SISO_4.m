import idModels.*

% Polynome A und B
A = conv([1 -.5],1); B = [0 .5]; y0 = 10; p1 = 15;

na = length(A)-1; nk = find(B~=0,1,'first')-1; nb = length(B)-nk; 

% % % Nichtlinearität f(u,y) = u1^2y+p1u1^3 definieren
F = @(u,p,y) deal(  u.^2.*y + u.^3.*p(1),...
                    u.^3,...
                    2*y.*u + 3*p(1).*u.^2,...
                    u.^2);
I.parameters = 15;
I.input_idx = 1;
I.output_idx = 1;
I.fun = F;

% Modell simulieren
u1 = kron(rand(10,1),ones(10,1));
%p=[p1,p2];
x(1,1) = 0; y = y0;
for k = 2:length(u1)
     [fu,~,~,~] = F(u1(k-1:-1:k-nb),p1,y(k-1:-1:k-nb));
     x(k,1) = -A(2:end)*x(k-1:-1:k-na) + B(2:end)*fu;
     y(k,1) = y0 + x(k,1);
end

% Modell identifizieren
M = idModels.NsfPolyModel(1,1,1,'InputNonlinearity',I);
M.identify(y,u1,'EstimateOutputOffset',1); 

% Modell Simulieren
[ys,x,dy_du] = M.simulate(u1);    %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
plot(y); hold on; plot(ys,'--r')

%% Check with finite differences
figure('color','w'); 
j = 10; %Zeitschritt für u für welchen der gradient berechnet werden soll!
delta = 1e-8;
for k = 1:size(u1,2)
    up = u1; um = u1;
    up(j,k) = up(j,k) + delta; 
    um(j,k) = um(j,k) - delta; 
    grad = (M.simulate(up) - M.simulate(um))/(2*delta); % numerical gradient
    Dif = grad - squeeze(dy_du(:,1,k,j)); % should be very small (>1e-6)
    subplot(size(u1,2),1,k); 
    plot(Dif);
    %plot(squeeze(dy_du(1,k,:,j)));hold on; plot(grad,'--r')
    ylabel(['dy/du(' num2str(j) ')'])
end
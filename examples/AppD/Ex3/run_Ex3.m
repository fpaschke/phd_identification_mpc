close all; clc; clear;
addpath('..\..\..');
import idModels.*

%% Systemmodell definieren und simulieren
% Struktur für f(u,y,eta) definieren
F(1).fun = 'f_1';                                 % Funktionsname von f_1
F(1).parameters = 2;                              % eta_1
F(1).input_idx = [1 2];                           % u1 und u2 an f_1
F(1).output_idx = 1;                              % y1 and f_1 übergeben
F(2).fun = 'f_2';                                 % Funktionsname von f_3
F(2).parameters = 5;                              % eta_2
F(2).input_idx = 2;                               % u2 an f_2 übergeben 
F(2).output_idx = 2;                              % y2 and f_2 übergeben

% System S instanziieren und Polynomparameter setzen
S = NsfPolyModel(ones(2),eye(2),eye(2),...        % System instanziieren
    'InputNonlinearity',F,'NoiseVariance',zeros(2));
S.A(1,1).val = [1 -.9]; S.A(1,2).val = [0 .5];    % A Parameter setzen 
S.A(2,1).val = [0 .1]; S.A(2,2).val = [1 -.5];
S.B(1,1).val = [0 .1]; S.B(2,2).val = [0 .05];    % B Parameter setzen

% System S Simulieren
rng(0);                                           % für Reproduzierbarkeit     
u = kron(rand(10,2),ones(50,1));                  % Inputs erzeugen
x = S.simulate(u);                                % Simulieren 
y = x + [.05 .02].*randn(size(u));                % Ausgangsrauschen add.

%% Modell definieren
% Modell M instanziieren, identifizieren und Parameter ausgeben
F(1).parameters = 1; F(2).parameters = 1;         % Startwerte Optimierung
M = NsfPolyModel(ones(2),eye(2),...               % Modell M instanziieren
    eye(2),ones(2,1),'InputNonlinearity',F);
M.identify(y,u);                                  % M identifizieren
M.printParameters;                                % Parameter ausgeben                        
[ys,~,dy_du] = M.simulate(u);                     % M simulieren

%% Plot
figure('Color','w','Position',[200 200 800 450]);    
s(1)=subplot(2,1,1); plot([y(:,1) ys(:,1)]); grid on; ylabel('$y_1$'); legend({'$y_1[t]$' '$\hat{y}_1[t]$'},'Location','southeast'); 
yyaxis right; L = plot(u(:,1)); ylabel('$u_1$'); s(1).XTickLabel = {}; L.Annotation.LegendInformation.IconDisplayStyle = 'off';
s(2)=subplot(2,1,2); plot([y(:,2) ys(:,2)]); grid on; ylabel('$y_2$'); legend({'$y_2[t]$' '$\hat{y}_2[t]$'},'Location','southeast');
yyaxis right; L = plot(u(:,2)); ylabel('$u_2$'); xlabel('$t$ in Zeitschritten'); L.Annotation.LegendInformation.IconDisplayStyle = 'off';
util.formatFigure(14);
util.saveTightFigure(gcf,'D:\Diss\Bilder\Anhang\Ex3.pdf','AxPosOffset',[0 0.09 0 0],'SubPlotYspace',.08,'FigPosOffset',[0 0 0 -40],'EqualX',true)

%% Gradient mit finiten differenzen überprüfen
figure('color','w'); 
j = 50; % 
delta = 1e-8;
m=1;
for l=1:size(y,2)
    i=1;
    for k = 1:size(u,2)
        up = u; um = u;
        up(j,k) = up(j,k) + delta; 
        um(j,k) = um(j,k) - delta;
        grad = (M.simulate(up) - M.simulate(um))/(2*delta); % numerical gradient
        Dif = grad - squeeze(dy_du(:,:,i,j)); % should be very small (>1e-6)
        subplot(size(u,2),size(y,2),m); 
        plot(Dif(:,l)); hold on; ylabel(['\Delta dy_{' num2str(l) '}/du_{' num2str(k) '}(' num2str(j) ')']);
        yyaxis right; plot(squeeze(dy_du(:,l,i,j))); ylabel(['dy_{' num2str(l) '}/du_{' num2str(k) '}(' num2str(j) ')']);
        m = m + 1;
        i = i + 1;  
    end   
end
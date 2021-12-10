clc; clear; close all;
addpath('..\..\..');
import idModel.*;

%% System deklarieren und simulieren
% Parameter und Signale definieren
rng(0);                                      % für Reproduzierbarkeit
T_d = 50;                                    % Periodendauer d_1
N = 500;                                     % Anz. Samples
u = 2*kron(rand(N/T_d/2,1),...               % Eingangssignal u
	[zeros(T_d,1); ones(T_d,1)]);
d1 = 1*kron(ones(N/T_d,1), ...               % periodische Störung d_1
     [zeros(T_d/2,1);ones(T_d/2,1)]);
d2 = cos(2*pi*linspace(0,1,N)')+.5;          % Störung d_2 ("Drift")    
A = [1 -.9]; nA = length(A)-1;               % Polynomkoff. A
B = [0 .5]; nB = length(B)-1;                % Polynomkoff. B
f = @(u,eta,y) deal(u.*(eta-y),...           % f(u,eta,y)
                    u);                      % df/deta(u,eta,y)
eta = 3;                                     % Parameter f(u,y,eta)     

% System Simulieren
x = NaN(N,1); y = x; x(1) = 0; y(1) = d2(1); % Initialisierung   
for t = 2:N
     [fuy,~] = f(u(t-1:-1:t-nB),...          % f(u,eta,y) auswerten
                   eta,y(t-1:-1:t-nB)); 
     x(t) = -A(2:end)*x(t-1:-1:t-nA) + ...   % Ausgang lin. Teilsystem
         B(2:end)*(fuy + d1(t-1:-1:t-nB));
     y(t) = x(t) + d2(t);                    % Systemausgang mit Drift      
end

%% Model deklarieren und identifizieren
% Struktur für Eingangsnichtlin. definieren
F(1).fun = f;                       % Funktionshandle f(u,eta,y) 
F(1).parameters = 0;                % Optimiererstartwert für eta
F(1).free = 1;                      % als freien Param. festlegen
F(1).input_idx = 1;                 % Ein- bzw. Ausgänge die bei 
F(1).output_idx = 1;                % Berechnung von f benötigt werden 

% Modell definieren, identifizieren und Parameter ausgeben
M = idModels.NsfPolyModel(...       % ARARX Modellansatz M definieren
      nA,1,nB,0,T_d,'InputNonlinearity',F);
M.D.val(2:end-1) = 0;               % d_1 bis d_(T_d-1) Null setzen 
M.D.free(2:end-1) = false;          % und als feste Param. festlegen
M.identify(y,u,'IntegrateNoise',1); % M mit Vorfilter L=1-q^-1 ident.
M.printParameters;                  % Parameter ausgeben

%% Plot
k = 30;
yp = M.calcPredictions(y,u,'Hp',k,'Nmax',length(M.Ss.A));
h = figure('Color','w','Position',[200 200 700 600]); col = lines(3);
s(1) = subplot(2,1,1); plot(y,'LineWidth',1); hold on; plot(x,'LineStyle',':','Color',col(1,:),'LineWidth',1); hold on; plot(squeeze(yp(:,:,end)),'LineWidth',1,'LineStyle','--','Color',col(2,:)); s(1).XTickLabel = {}; ylim([-.5 6.5]);
s(2) = subplot(2,1,2); plot(u,'LineWidth',1); hold on; plot(d1,'LineStyle','--','Color',col(2,:),'LineWidth',1); hold on; plot(d2,'LineStyle',':','Color',col(2,:),'LineWidth',1); 
util.formatFigure(14,{'' '$t \textrm{ in Zeitschritten}$'},{'Ausgang' 'Eingang/St\"orung'},[],[],{{'$y[t]$' '$\tilde{y}[t]$' ['$\hat{y}[t|t-' num2str(k) ']$']} {'$u[t]$' '$v_1[t]$' '$v_2[t]$'}},'LegLoc','SouthEast')
%util.saveTightFigure(h,'D:\Diss\Bilder\Identifikation\BilSys_sim.pdf','SubPlotYspace',.08,'FigPosOffset',[0 0 0 -80],'AxPosOffset',[0.0 0.06 0 0],'EqualX',1);
clc; clear; close all;
addpath('..\..\..');

%% Generate Data and simulate System
rng(0)                                    % für Reproduzierbarkeit                            
u = kron(rand(20,1),ones(50,1));          % generiere Eingangssignal
B = [0 1e-1];                             % Zählerkoeffizienten
A = conv([1 -.95],[1 -.85]);              % Nennerkoeffizienten
x = filter(B,A,u) + 3;                    % simuliere System
y = round(x);                             % quantisierter Ausgang

M = idModels.NsfPolyModel(2,1,1,2);       % ARMAX model erzeugen
M.factorize('A');                         % faktorisiere A(q^-1) 
M.identify(y,u,'EstimateOutputOffset',1); % PEM Identifkation
M.printParameters;                        % Param. auf Konsole ausgeben
M.show('Prediction',y,u,'Hp',100);        % Plot Messung vs. Prädiktion

%% Output
f = gcf;
subplot(2,1,1); ylabel('$y$'); legend({'$y[t]$' '$\hat{y}[t|t-k]$'},'Location','southeast');
subplot(2,1,2); ylabel('$u$'); legend off; xlabel('Zeit $t$ in Zeitschritten','Interpreter','latex')
util.formatFigure(13);
util.saveTightFigure(f,'D:\Diss\Bilder\Anhang\Ex1.pdf','AxPosOffset',[0 0.09 0 0],'SubPlotYspace',.08,'FigPosOffset',[0 0 0 -40])

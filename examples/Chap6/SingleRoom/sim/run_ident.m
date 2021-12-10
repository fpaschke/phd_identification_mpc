% Generates data for identification by simulation,
% estimates models and saves them into resource folder.

clear; clc;               
addpath('..\..\..\..');
addpath('..\model\')
addpath('..\..\LectureBuilding\sim\subfun\sig\')
addpath('..\..\LectureBuilding\sim\subfun\ident')
addpath(genpath('.\subfun\'))

%% Estimate Room Model
Hp_val_flo = 12;                             	% Time Horizon for validation [h]
P = getParameters('T_s',15,'Hp_val',Hp_val_flo,'PeriodicOccupancy',true);    % Generate default Parameters
P.Sig.N_pers.NumPeop = 0;                       % Use some number if seasonal model should be estimated

P.Tstart = '2020-3-1';                          % Starting Time of Simulation
P.N_w = 8;                                   	% Number of simulation [weeks] (length of training data)
SimOutFlo = simulate('SingleRoomOpenLoop',P); 	% Simulate Open Loop

% Identify Heating Curve
HcvModFlo = estHcv(SimOutFlo);

% Identify Model with Prefiltering (use P.Ident to set identification options)
P.Ident.fc = 1/(2*6*pi);                        % Cutoff Frequency of highpass butterworth prefilter in 1/h. 
P.Ident.Seasonality = 672;                    	% Use 672 for weekly seasonality                          
P.Ident.Orders = {2 [2 1 1] 1 2 0 0};       	% Model orders na,nb,nk,nc,nd,ne
RoomModFlo = estRoom(SimOutFlo,P.Ident);        % Estimate Room Model
%ylim([.2 .3]); util.saveTightFigure(gcf,'D:\Diss\Bilder\MPC\SR_RMMSE.pdf','AxPosOffset',[0 0 -.01 0])

%% Estimate S-ARMA Weather Models using Multistep Criterion (Ts = 1h)
P.Ident.Orders = 1:3;                           % Order of Polys = 1
P.Ident.HpVal = Hp_val_flo;                   	% Horizon for validation
P.Ident.fc = [];                                % Cutoff Frequency of highpass butterworth prefilter in 1/h. 
P.Ident.FactorizeG = false;                     % Do not factorize A
P.Ident.Seasonality = 24;                       % Seasonality of Model
P.Ident.SeasonalityOrder = 'ByOrder';  			% Order of seasonal part (pos int.), if 'ByOrder' then it will be increased with P.Ident.Orders 
P.Ident.Iopt.ForcePolesG = '';                  % No restrictions on poles of G
P.Ident.Iopt.IntegrateNoise = false;            % No noise integration
P.Ident.Iopt.ForceSymmetricA = false;           % Do not force symmetric A
P.Ident.Iopt.Algorithm = 'levenberg-marquardt'; % Optimizer alg
WeaMod = estWeatherModel(In,P.Ident);           % Estimate Weather Model. 

%% Show and validate Roommodel
RoomModFlo.TimeUnit = 'hours'; RoomModFlo.Ts = 0.25;
RoomModFlo.showBode(); set(gcf,'Position',[100 100 1200 300])                       % Plot Bode Diagramm    
subplot(1,4,1); ylabel(''); for i = 1:3; set(subplot(1,4,i),'XTick',[1e-4 1e-2 1]); end
s = subplot(1,4,4); cla reset; RoomModFlo.showInputNonlinearity(3,[NaN 1],0,'Axes',s,'Color','k'); xlabel('Ventilpos. $H$ [\%]'); ylabel(''); title('Ventilkennlinie $H^{\eta_1}$','Interpreter','latex'); 
util.formatFigure(14); set(gca,'XTickLabel',num2cell(100*cellfun(@str2num,get(gca,'XTickLabel'))))
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\MPC\SR_Bode.pdf','XShift',-.05,'FigPosOffset',[0 0 -130 30],'AxPosOffset',[-0.07 0.065 0 0],'SubplotXSpace',-0.02)
RoomModFlo.TimeUnit = 'seconds'; RoomModFlo.Ts = P.T_s*60;

% Get data
Y = SimOutFlo.T_air.Data; N = size(Y,1);
U =  [SimOutFlo.In.T_out.Data(1:N) ...
      SimOutFlo.In.P_sun.Data(1:N) ...
      SimOutFlo.In.Azimuth.Data(1:N) ...
      SimOutFlo.In.Elevation.Data(1:N) ...
      SimOutFlo.H_win.Data(1:N) ...
      SimOutFlo.H_flo.Data(1:N) ...
      SimOutFlo.T_sup_flo.Data(1:N)];
  
figure('color','w','Position',[100 100 750 270]);
s(1) = subplot(1,2,1); RoomModFlo.show('Boxplot',Y,U,'Hp',Hp_val_flo*4,'Axes',s(1),'DeltaK',2); legend('off')
s(2) = subplot(1,2,2); RoomModFlo.show('Autocovariance',Y,U,'MaxLag',4*30,'Significance',.01,'Axes',s(2)); legend('off')
s(1).XAxis.TickValues = 0:4:24; s(1).XAxis.TickLabels = {'0' '2' '4' '6' '8' '10' '12'}; s(1).XLabel.String = 'Pr\"adiktionshorizont [h]'; s(1).YLabel.String = '$\varepsilon[t]$';
s(2).XAxis.TickValues = 0:24:120; s(2).XAxis.TickLabels = {'0' '6' '12' '18' '24' '30'}; s(2).XLabel.String = 'Lag $\tau$ [h]'; s(2).YLabel.String = 'Autokovarianz $\hat{r}_\varepsilon[\tau]$';
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\MPC\SR_BoxAcf.pdf','XShift',-.05,'FigPosOffset',[0 0 -75 20],'AxPosOffset',[0 .06 0 0])


Pv = SimOutFlo.H_flo.Data(1:N).^RoomModFlo.InputNonlinearity(3).parameters{1}.*(SimOutFlo.T_sup_flo.Data(1:N) - SimOutFlo.T_air.Data(1:N)).^1.1;
figure('color','w'); scatter(SimOutFlo.P_flo.Data,Pv,12,'filled','MarkerEdgeAlpha',.1,'MarkerFaceAlpha',.1); grid on; xlabel('Heizleistung'); ylabel('Virtuelle Heizleistung');

%% Plot Models
% Generate Predictions
ypred = RoomModFlo.calcPredictions(SimOutFlo.T_air.Data,U,'Hp',Hp_val_flo*3600/RoomModFlo.Ts,'ObserverInitSamples','auto'); 
Tpred = SimOutFlo.T_air; 

figure('color','w','Position',[100 100 900 750]); 
s(1) = subplot(3,1,1); 
plot(SimOutFlo.T_air,'DisplayName','Modelica Model'); hold on;
%Tpred.Data = ypred(:,:,end/4); plot(Tpred,'DisplayName',[num2str(Hp_val_flo/4) 'h Pr\"adiktion']); hold on; 
Tpred.Data = ypred(:,:,end/2); plot(Tpred,'DisplayName',[num2str(Hp_val_flo/2) 'h Pr\"adiktion']); hold on; 
%Tpred.Data = ypred(:,:,end); plot(Tpred,'DisplayName',[num2str(Hp_val_flo) 'h Pr\"adiktion']); hold on; 
ylabel('$\vartheta_\mathrm{R}$ [$^{\circ}$C]'); ylim([19.0 26]); legend('Location','NorthWest'); grid on; title(''); set(gca,'XTickLabel',{})
s(2) = subplot(3,1,2); 
plot(100*SimOutFlo.H_flo,'DisplayName','Heizventil'); hold on; ylim([-5 105]); ylabel('$H$ [\%]'); yyaxis right;
plot(100*SimOutFlo.H_win,'DisplayName','Verschattung'); hold on; ylim([-5 105]); ylabel('$H_\mathrm{F}$ [\%]');
grid on; title(''); set(gca,'XTickLabel',{}); 
s(3) = subplot(3,1,3);
plot(SimOutFlo.In.T_out); ylabel('$\vartheta_\mathrm{Aul}$ [$^{\circ}$C]'); grid on; yyaxis right;
plot(SimOutFlo.In.P_sun/1e3); ylabel('$\dot{Q}_\mathrm{GS}$ [kW/m$^2$]'); grid on; ylim([-0.05 0.6]);
title(''); linkaxes(s,'x'); xlim(datetime({'2020-3-1' '2020-3-28'})); util.formatFigure(14);
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\MPC\SR_Pred.pdf','XShift',-.05,'FigPosOffset',[0 0 10 -160],'AxPosOffset',[0.055 -.015 0 0],'SubplotYspace',0.05,'EqualX',true)

%% Plot Weather Model
figure('color','w'); 
s(1) = subplot(2,1,1);
plot(In.T_out,'-k','DisplayName','Außentemp. [C]'); hold on; xl = xlim;
ypred = ToutMod.calcPredictions(In.WeatherData.data(:,2),[],Hp_val_flo*3600/ToutMod.Ts,[],[],[]);
Pred = timeseries(ypred(:,1,end),0:(length(ypred)-1)); Pred.TimeInfo.StartDate = '2020-1-1'; Pred.TimeInfo.Units = 'hours';
plot(Pred,'--r','DisplayName',['Prädiktion Außentemp. [C] (' num2str(P.Ident.HpVal) '-Schritt)']); grid on; xlim(xl);
ylabel('Außentemp. [C]'); title('');
s(2) = subplot(2,1,2);
plot(In.P_sun,'-k','DisplayName','Globalstrahlung. [kW/m^2]'); xlim(xl); grid on
ypred = PsunMod.calcPredictions(In.WeatherData.data(:,9)/1e3,[],Hp_val_flo*3600/PsunMod.Ts,[],[],[]); hold on;
Pred.Data = ypred(:,1,end);
plot(Pred,'--r','DisplayName',['Prädiktion Globalstrahlung. [kW/m^2] (' num2str(P.Ident.HpVal) '-Schritt)']); xlim(xl); grid on;
ylabel('Globalstrahlung. [kW/m^2]'); title('');
linkaxes(s,'x');

%% Save Models
save(['.\resource\S-Armax_Ts'  num2str(P.T_s) '_Hp15m'],...
    'RoomModFlo','HcvModFlo','WeaMod');
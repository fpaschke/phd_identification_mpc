% Generates data for identification by simulation of open loop,
% estimates models and saves them into resource folder.

clear; clc;        
addpath('..\..\..\..');
addpath('..\model\')
addpath('..\..\SingleRoom\sim\subfun\sig\')
addpath(genpath('.\subfun\'))

%% Estimate Room Models with Seasonality
Hp_val = 3*24;                                  	% Time Horizon for validation [h]
P = getParameters('T_s',5,'T_s_ident',120,'Hp_val',Hp_val,'N_w',24);% Generate default Parameters
P.Sig.N_pers.Mean(:) = 0;                        	% No Occupancy

P.Tstart = '2021-1-1';                              % Starting Time of Simulation
SimOut = simulate('LecBuilOpenLoop',P);             % Simulate

%% Identify (use P.Ident to set identification options)
% Estimate 1 MISO - Model for each zone 
P.Ident.MIMO = false;                        	% Estimate MISO Models and merge them into one
P.Ident.Seasonality = 7*24/2;               	% Seasonality order (1 week = 7*24 Ts_ident  = 60): Use [] for no seasonality
P.Ident.SeasonalityOrder = 2;                   % Size of Seasonality noise model
P.Ident.fc = 1/(2*7*24*pi);                   	% Cutoff Frequency of highpass butterworth prefilter in 1/h.                 
P.Ident.Orders = {2 [2 2 2] 1 1};            	% nA, nBi nKi nC
RoomMod_SO_Seas = estRoomModel(SimOut,P.Ident);      % Estimate Room Models

% Estimate coupled MIMO - Model (Non Factorized polynomials)
P.Ident.MIMO = true;                            % Estimate MIMO Model
P.Ident.Seasonality = [];                   	% Seasonality order (1 week) 
P.Ident.SeasonalityOrder = 2;                   % Size of Seasonality noise model
P.Ident.Iopt.ForceSymmetricA = true;         	% Force symmetric A Matrix of ARMAX Model 
P.Ident.fc = [];                             	% No Prefiltering
P.Ident.Iopt.Algorithm = 'sqp';
nA =   [2 0 1 0 1; 
        0 2 0 1 1; 
        1 0 2 0 1;
        0 1 0 2 1;
        1 1 1 1 2];
nB =   [    2 0 0 0 2 0 0 0 0 0 2; 
            0 2 0 0 0 2 0 0 0 0 2; 
            0 0 2 0 0 0 2 0 0 0 2;
            0 0 0 2 0 0 0 2 0 0 2;
            0 0 0 0 0 0 0 0 2 2 2];    
nC = 1*ones(5,1);
P.Ident.Orders = {nA nB ones(size(nB)) nC};
RoomMod_MO = estRoomModel(SimOut,P.Ident); % Estimate MIMO Model in ARMAX Structure

%% Estimate MO-S-ARMA Weather Model using 1-Step PEM 
P.Ident.Orders = 1:3;                           % Order of Polys = 1
P.Ident.HpVal = Hp_val;                         % Horizon for validation
P.Ident.fc = [];                                % Cutoff Frequency of highpass butterworth prefilter in 1/h. 
P.Ident.FactorizeG = false;                     % Do not factorize A
P.Ident.Seasonality = 24;                       % Seasonality of Model
P.Ident.SeasonalityOrder = 'ByOrder';  			% Increase Order of Seasonal part
P.Ident.Iopt.IntegrateNoise = false;
P.Ident.Iopt.ForceSymmetricA = false;
P.Ident.Iopt.Algorithm = 'levenberg-marquardt';
WeaMod = estWeatherModel(In,P.Ident);           % Estimate Weather Model.        

%% Save Models
save(['.\resource\Models_Ts' num2str(P.Ident.T_s)],'RoomMod_MO','RoomMod_SO','RoomMod_SO_Seas','WeaMod');

%% Plot Room Models
% Generate Predictions
[N,Nz] = size(SimOut.T_air.Data);
dk = RoomMod_SO.Ts/SimOut.T_air.TimeInfo.Increment;
U =  [SimOut.In.T_out.Data(1:dk:N) ...
      SimOut.In.P_sun.Data(1:dk:N) ...
      SimOut.In.Azimuth.Data(1:dk:N) ...
      SimOut.In.Elevation.Data(1:dk:N) ...
      SimOut.T_sup.Data(1:dk:N) ...
      SimOut.m_sup.Data(1:dk:N,:)];
Y = SimOut.T_air.Data(1:dk:N,:);  
ypred_SO = RoomMod_SO.calcPredictions(SimOut.T_air.Data(1:dk:N,1:Nz-1),U,'Hp',Hp_val*3600/RoomMod_SO.Ts); 
Pred_SO = SimOut.T_sup.resample(SimOut.T_sup.Time(1:dk:N));
Pred_SO.Data = ypred_SO(:,:,end);
ypred_MO = RoomMod_MO.calcPredictions(SimOut.T_air.Data(1:dk:N,:),U,'Hp',Hp_val*3600/RoomMod_SO.Ts);
Pred_MO = SimOut.T_sup.resample(SimOut.T_sup.Time(1:dk:N));
Pred_MO.Data = ypred_MO(:,:,end);

figure('color','w','Position',[100 100 900 1000]); C = lines(4); 
for nz = 1:Nz-1  
    S(nz) = subplot(Nz-1,1,nz); 
 	plot(getSubSignal(SimOut.T_air,nz),'-k','DisplayName','Modelica Modell'); hold on;
    plot(getSubSignal(Pred_SO,nz),'Color',C(1,:),'DisplayName',[num2str(size(ypred_SO,3)*RoomMod_SO.Ts/3600) 'h Pr\"ad. (MISO)']); hold on;
    plot(getSubSignal(Pred_MO,nz),'Color',C(2,:),'DisplayName',[num2str(size(ypred_SO,3)*RoomMod_MO.Ts/3600) 'h Pr\"ad. (MIMO)']); hold on;
    plot(SimOut.T_sup,'Color',C(3,:),'DisplayName','Zuluft'); hold on;   
    ylabel(['$\vartheta_' num2str(nz) '$[$^{\circ}$C]']); grid on; title(''); set(gca,'XTickLabel',{}); ylim([16.0 25.5]);
    yyaxis right; 
    plot(getSubSignal(SimOut.m_sup,nz),'Color',C(4,:),'DisplayName','Massenstrom');
    if nz <= 2; ylim([-.5 6.5]); else; ylim([-.5 4.5]); end
    S(nz).YAxis(2).Color = C(4,:);
    ylabel(['$\dot{m}_' num2str(nz) '$ [kg$/$s]']); if nz ~= Nz-1; set(gca,'XTickLabel',{}); end 
end
linkaxes(S,'x'); util.formatFigure(14); 
%xlim(datetime({'2021-4-1' '2021-4-21'}));
xlim(datetime({'2021-2-27' '2021-4-1'}));
util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\MPC\LB_Pred_2h.pdf','XShift',-.05,'FigPosOffset',[0 0 10 -260],'AxPosOffset',[0.055 -.04 0 0],'SubplotYspace',0.04,'EqualX',true)

%% ACF/BOX-PLOT
figure('color','w','Position',[100 100 1100 450]);
for nz = 1:Nz-1 
    s(1,nz) = subplot(2,Nz-1,nz); 
    %dk = 3;
    dk = 2;
    RoomMod_SO.show('Boxplot',Y(:,1:Nz-1),U,'Hp',Hp_val*3600/RoomMod_SO.Ts,'Axes',s(1,nz),'DeltaK',dk,'Color',C(1,:),'IxOut',nz); legend('off'); hold on;
    RoomMod_MO.show('Boxplot',Y,U,'Hp',Hp_val*3600/RoomMod_SO.Ts,'Axes',s(1,nz),'DeltaK',dk,'Color',C(2,:),'IxOut',nz); legend('off');
    s(1,nz).XLabel.String = '$k$'; 
    %ylim([-.5 .5]);
    ylim([-.75 .75]);
    s(2,nz) = subplot(2,Nz-1,nz+Nz-1);
    RoomMod_SO.show('Autocovariance',Y(:,1:Nz-1),U,'MaxLag',8*24*3600/RoomMod_SO.Ts,'Significance',.01,'Axes',s(2,nz),'IxOut',nz,'Color',C(1,:)); legend('off');
    RoomMod_MO.show('Autocovariance',Y,U,'MaxLag',8*24*3600/RoomMod_SO.Ts,'Significance',.01,'Axes',s(2,nz),'IxOut',nz,'Color',C(2,:)); legend('off');
    %ylim([-2 10]*1e-3);
    ylim([-5 20]*1e-3);
    %s(2,nz).XTick = 0:24:72; s(2,nz).XTickLabel = {'0' '24' '48' '72'};
    s(2,nz).XTick = 0:24:8*12; s(2,nz).XTickLabel = {'0' '24' '48' '72' '96'};
    if nz == 1
      	s(1,nz).YLabel.String = ['$\varepsilon_i[t|t-k]$']; 
        s(2,nz).YLabel.String = ['$\hat{r}_{\varepsilon_i}[\tau]$'];
    else
        s(1,nz).YLabel.String = ''; 
        s(2,nz).YLabel.String = ''; 
        s(1,nz).YTickLabel = {};
        s(2,nz).YTickLabel = {};
    end
end
util.formatFigure(14);
for nz = 1:Nz-1
    s(1,nz).Position(1) = s(1,nz).Position(1) - .07 - (nz-1)*0.02; s(1,nz).Units = 'inches'; 
    s(2,nz).Position(1) = s(2,nz).Position(1) - .07 - (nz-1)*0.02; s(2,nz).Units = 'inches';
end
set(gcf,'Position',[100 100 880 430],'Units','Inches')
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'\\tsclient\D\Diss\Bilder\MPC\LB_BoxAcf_2h.pdf','-dpdf')

%% Plot Weather Models
figure('color','w'); 
s(1) = subplot(2,1,1);
ToutDat = In.WeatherData.data(:,2);
Tout = timeseries(ToutDat,0:(length(ToutDat)-1)); Tout.TimeInfo.StartDate = '2020-1-1'; Tout.TimeInfo.Units = 'hours';
plot(Tout,'-k','DisplayName','Außentemp. [C]'); hold on
ypred = ToutMod.calcPredictions(In.WeatherData.data(:,2),[],Hp_val*3600/ToutMod.Ts,[],[],[]);
Pred = timeseries(ypred(:,1,end),0:(length(ypred)-1)); Pred.TimeInfo.StartDate = '2020-1-1'; Pred.TimeInfo.Units = 'hours';
plot(Pred,'--r','DisplayName',['Prädiktion Außentemp. [C] (' num2str(P.Ident.HpVal) '-Schritt)']); grid on;
ylabel('Außentemp. [C]'); title('');
s(2) = subplot(2,1,2);
Psun = timeseries(In.WeatherData.data(:,9)/1e3,0:(length(ToutDat)-1)); Psun.TimeInfo.StartDate = '2020-1-1'; Psun.TimeInfo.Units = 'hours';
plot(Psun,'-k','DisplayName','Globalstrahlung. [kW/m^2]'); grid on
ypred = PsunMod.calcPredictions(In.WeatherData.data(:,9)/1e3,[],Hp_val*3600/PsunMod.Ts,[],[],[]); hold on;
Pred.Data = ypred(:,1,end);
plot(Pred,'--r','DisplayName',['Prädiktion Globalstrahlung. [kW/m^2] (' num2str(P.Ident.HpVal) '-Schritt)']); grid on;
ylabel('Globalstrahlung. [kW/m^2]'); title('');
linkaxes(s,'x');

clc; clear; close all;
addpath('..\..\..')

%% Load Data
params.Fname = 'Data_ConferenceRoom.xls'; % Name of data file
[Info, Sheets] = xlsfinfo(params.Fname); % Load xlsinfo to get number of sheets
D = cellfun(@(sh) xlsread(params.Fname,sh),Sheets(2:end)','UniformOutput',false); % Load Data to tables
[~,ID] = xlsread(params.Fname,'Data1','B1:X1'); % Load column names        

%% Set Id Options
params.Ts = xlsread(params.Fname,'DataInfo','B3');      % Sampling Time of Model 
params.Win = xlsread(params.Fname,'DataInfo','B6:B7');  % Orientation of Windows (wrt North and to a hor. plane)
params.HpVal = 1*60/params.Ts;                         	% Prediction horizon in #steps for validation
params.Orders = 1:6;                                 	% Orders of Polynomials A,B,... for which identification will be performed
params.Types = {'arx' 'armax' 'bj' 'oe'};            	% Parameterizations for which identification will be performed
params.InitFromPreviosN = true;                     	% Initialize Model of Order n with estimated Model with Order n-1
params.FactorizeG = true;                               % Identifiy model with factorized A(q),B(q),E(q)
params.fc = 1/(2*2*pi);                                 % Cutoff Frequency of highpass butterworth prefilter in 1/h. 
params.HpEst = [1 params.HpVal];                        % Prediction Horizon for estimation

% Options for Identification
params.iopt.TolX = 1e-8;                             	% Stops if norm of change in parametervector is less then TolX
params.iopt.ForcePolesG = '';                           % Forces stable pos Real Poles to G 
params.iopt.Hp = 1;                                     % Prediction horizon in #steps for identification
params.iopt.MaxFunEval = 1e3;                           % Max number of Cost function calls.  
params.iopt.MaxIter = Inf;                              % Max number of iterations during optimization.  
params.iopt.FunTolAbs = 1e-8;                       	% Abort criterion of optimization. Abort if |F_i - F_i-1| <= FunTolAbs
params.iopt.StabilizationMethod = 'none';               % Stabilization of opt. predictor
params.iopt.Ic = 'backcast';                            % Initialization Method of residual filter 
params.iopt.IntegrateNoise = false;                   	% Force intergrator in noise model. 
params.iopt.InitMethod = 'arx';                         % Initialize Model with least squares
params.iopt.Solver = 'matlab_lsqnonlin';                % Chosen optimization function. If no optimization toolbox installed you can use 'opti_levmar' from OPTI-Toolbox instead. (See https://inverseproblem.co.nz/OPTI/)
params.iopt.Algorithm = 'levenberg-marquardt';          % Chosen optimization algorithm (Remove this line if Solver = opti_levmar)
params.iopt.MultiStep = true;                       	% Use Multistep Criterion for identification  if Hp>1 (see below)    

%% HAMMERSTEIN MODEL: Extract Data and Identify
% Model Setup
params.col_in = [2 3 4 5 10 11]; % Outsidetemp., global radiation on horizontal plane, sun height, sun azimuth, Heating power, cooling power
params.col_out = 1; % Roomtemperture 
params.mdl = struct;
params.mdl.Name = 'Hammerstein-LTI';
params.mdl.Ts = params.Ts;
params.mdl.TimeUnit = 'minutes';
params.mdl.InputName = ID(params.col_in); 
params.mdl.OutputName  = ID(params.col_out);
params.mdl.InputNonlinearity(1).fun = ''; % Identity
params.mdl.InputNonlinearity(1).parameters = [];
params.mdl.InputNonlinearity(1).free = []; 
params.mdl.InputNonlinearity(1).input_idx = 5; % Heating Power
params.mdl.InputNonlinearity(2).fun = ''; % Identity
params.mdl.InputNonlinearity(2).parameters = [];
params.mdl.InputNonlinearity(2).free = []; 
params.mdl.InputNonlinearity(2).input_idx = 6; % Cooling Power
params.mdl.InputNonlinearity(3).fun = ''; % Identity
params.mdl.InputNonlinearity(3).parameters = [];
params.mdl.InputNonlinearity(3).free = []; 
params.mdl.InputNonlinearity(3).input_idx = 1; % OutsideTemp
params.mdl.InputNonlinearity(4).fun = 'idModels.func.fun_rad'; % calculates radiation on inclined plane and multiply by logistic function to consider shading by trees and buildings in front of the window
params.mdl.InputNonlinearity(4).parameters = [params.Win' 2]; % orientation of windows wrt N in degrees, orientation of windows wrt to horizontal plane in degrees, method for computing diffuse radiation on inclined plane. see util.calcRadiationOnInclinedPlane,estimated shading angle, growth rate/steepness of logistic function
params.mdl.InputNonlinearity(4).free = [false false false]; 
params.mdl.InputNonlinearity(4).input_idx = [2 4 3]; % Radiation, Azimuth, Height   
params.mdl.PreFilter = cell(1,2); 
[params.mdl.PreFilter{:}] = butter(1,(params.fc/60)/((1/params.Ts)/2),'high');
params.mdl.PreFilter{1} = params.mdl.PreFilter{1}./params.mdl.PreFilter{1}(1);

% Extract Data
[y,u1] = cellfun(@(S) deal(S(:,params.col_out),S(:,params.col_in)),D,'UniformOutput',false); 

% Do identification for different Prediction horizons and Model Structures and plot (Abb. 5.2)
for k = 1:length(params.HpEst)
    params.iopt.Hp = params.HpEst(k);                     
    [RES_hamm{k}, M_hamm{k}] = idModels.test.sweepPar({'Order' 'Type'},@idModels.NsfPolyModel.estModel,y,u1,'Order',params.Orders,...
                                        'Type',params.Types,'InitFromPreviousN',params.InitFromPreviosN,'MdlProp',params.mdl,...
                                        'IdOpt',params.iopt,'HpVal',params.HpVal,'FactorizeG',params.FactorizeG,'MinRealTol',1e-6);
end

figure('color','w','Name','HAMMERSTEIN','Position',[100 100 1000 320]);
for k = 1:length(RES_hamm)
    [mat, lab] = util.tab2mat(RES_hamm{k},1); subplot(1,length(RES_hamm),k); bar(mat); grid on; ylim([0.14 0.26]);
    [mat, lab] = util.tab2mat(RES_hamm{k},1); subplot(1,length(RES_hamm),k); bar(mat); grid on; ylim([0.14 0.26]); 
end
util.formatFigure(14,'Polynomordnungen $n$',{'RMMSE($\varepsilon[t]$)' ''},[],[],{'ARMAX' 'ARX' 'BJ' 'OE'},'LegLoc','NorthEast');
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\Anwendungen\R51_Hamm_Cost.pdf','XShift',-.05,'FigPosOffset',[0 0 -90 25],'AxPosOffset',[0 .065 0 0]);

% Estimate ARIMAX Model with specified model orders 
params.iopt.MultiStep = false; 
nA = 2;
nB = [1 2 1 1];
nC = 1;
for k = 1:length(params.HpEst)
  	M_arimax{k} = idModels.NsfPolyModel(nA,nB,ones(size(nB)),nC,params.mdl);
  	M_arimax{k}.factorize({'A' 'B'});
    if k>1 % Init with previous Model
        M_arimax{k}.initParameters(M_arimax{k-1});
    end
    params.iopt.Hp = params.HpEst(k);
    M_arimax{k}.identify(y,u1,params.iopt);
    M_arimax{k}.updateSs('MinrealTol',1e-6); % Do PZ-Cancellation of StateSpace Object
	E = cell2mat(M_arimax{k}.calcResiduals(y,u1,params.HpVal));
	E = E(:); E = E(~isnan(E));
    RMMSE(k) = sqrt(mean(E.^2)); % RMMSE ERROR
end
T_1 = idModels.util.calcTimeconstant(M_arimax{1}.A.val(2:end),params.mdl.Ts)
T_12 = idModels.util.calcTimeconstant(M_arimax{2}.A.val(2:end),params.mdl.Ts)

% Plot Bode Diagramm 
M_arimax{1}.showBode('Models',M_arimax{2},'Legend',{},'ShowFrequencyWeighting',true,'ShowH',false);
set(gcf,'Position',[100 100 1200 600])
util.formatFigure(14);
for na = 1:5
    s = subplot(2,5,na);
    if na == 4; ylim([-120 -40]); elseif na == 2; ylim([-30 10]); end
    s.Position(2) = s.Position(2)-.08;
end    
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\Anwendungen\R51_Hamm_Bode.pdf','FigPosOffset',[0 0 -170 -70],'AxPosOffset',[-.07 -0.02 0 0]);

% Plot BoxPlot and ACF
h1 = figure('color','w','Name','BoxAcf','Position',[100 100 800 300]); 
s(1)=subplot(1,2,1); col = lines(2);
M_arimax{1}.show('Boxplot',y,u1,'Hp',params.HpEst(end),'Axes',s(1),'Color',col(1,:)); legend('off');  
M_arimax{2}.show('Boxplot',y,u1,'Hp',params.HpEst(end),'Axes',s(1),'Color',col(2,:)); legend('off');ylim([-.3 .3])
s(2)=subplot(1,2,2); 
M_arimax{1}.show('Autocovariance',y,u1,'MaxLag',72,'Significance',.01,'Axes',s(2),'color',col(1,:)); legend('off');
M_arimax{2}.show('Autocovariance',y,u1,'MaxLag',72,'Significance',.01,'Axes',s(2),'color',col(2,:)); legend('off')
util.formatFigure(14);
%util.saveTightFigure(h1,'\\tsclient\D\Diss\Bilder\Anwendungen\R51_Hamm_BoxAcf.pdf','AxPosOffset',[-0.04 .06 0 0],'FigPosOffset',[0 0 -60 25])

%% NONLINEAR MODEL (with supply temp. Measurement): Extract Data and Identify
% Model Setup
params.col_in = [2 3 4 5 7:9 12 6 13]; % Outsidetemp., global radiation on horizontal plane, sun height, sun azimuth, heating valve positions, heating supply temp, fan level, cooling supply temp.
params.col_out = 1; % Roomtemperture 
params.mdl = struct;
params.mdl.Name = 'NSF-LTI';
params.mdl.Ts =  params.Ts;
params.mdl.TimeUnit = 'minutes';
params.mdl.InputName = ID(params.col_in); 
params.mdl.OutputName  = ID(params.col_out);
params.mdl.InputNonlinearity(1).fun = 'idModels.func.fun_powdT'; % "virt. Heatingpower": mean(Valvepos)^p1  * (Tsupply_heat - Troom)^p2 
params.mdl.InputNonlinearity(1).parameters = [1 1.3]; % initial estimates of p1 and p2
params.mdl.InputNonlinearity(1).free = [true false]; 
params.mdl.InputNonlinearity(1).input_idx = 5:8; % Valve 1-3, SupplyTemp
params.mdl.InputNonlinearity(1).output_idx = 1; % Feedback of room temp
params.mdl.InputNonlinearity(2).fun = 'idModels.func.fun_powdT'; % "virt. Coolingpower": FanLvL^p1  * (Tsupply_cool - Troom)^p2 
params.mdl.InputNonlinearity(2).parameters = [1 1.0]; % initial estimates of p1 and p2
params.mdl.InputNonlinearity(2).free = [true false]; 
params.mdl.InputNonlinearity(2).input_idx = 9:10; % FanLvl, Cooling
params.mdl.InputNonlinearity(2).output_idx = 1; % Feedback of room temp
params.mdl.InputNonlinearity(3).fun = ''; % Identity function for outside temp.
params.mdl.InputNonlinearity(3).parameters = [];
params.mdl.InputNonlinearity(3).free = []; 
params.mdl.InputNonlinearity(3).input_idx = 1; % OutsideTemp
params.mdl.InputNonlinearity(4).fun = 'idModels.func.fun_rad'; % calculates radiation on inclined plane and multiply by logistic function to consider shading by trees and buildings in front of the window
params.mdl.InputNonlinearity(4).parameters = [params.Win' 2]; % orientation of windows wrt N in degrees, orientation of windows wrt to horizontal plane in degrees, method for computing diffuse radiation on inclined plane. see util.calcRadiationOnInclinedPlane,estimated shading angle, growth rate/steepness of logistic function
params.mdl.InputNonlinearity(4).free = [false false false]; 
params.mdl.InputNonlinearity(4).input_idx = [2 4 3]; % Radiation, Azimuth, Height     
params.mdl.PreFilter = cell(1,2); 
[params.mdl.PreFilter{:}] = butter(1,(params.fc/60)/((1/params.Ts)/2),'high');
params.mdl.PreFilter{1} = params.mdl.PreFilter{1}./params.mdl.PreFilter{1}(1);

% Extract Data
[~,u2] = cellfun(@(S) deal(S(1:params.Ts/5:end,params.col_out),S(1:params.Ts/5:end,params.col_in)),D,'UniformOutput',false); 

% Do identification for different Ic Methods
params.iopt.Hp = 1;	% Prediction horizon in #steps for identification
[RES_nl,M_nl] = idModels.test.sweepPar({'Order' 'Type'},@idModels.NsfPolyModel.estModel,y,u2,'Order',params.Orders,...
                                    'Type',params.Types,'InitFromPreviousN',params.InitFromPreviosN,'MdlProp',params.mdl,...
                                    'IdOpt',params.iopt,'HpVal',params.HpVal,'FactorizeG',params.FactorizeG,'MinRealTol',1e-6);

figure('color','w','Name','NSFLTI','Position',[100 100 500 330]);
[mat, lab] = util.tab2mat(RES_nl,1); bar(mat); ylim([0.16 0.24]); grid on;
util.formatFigure(15,'Polynomordnungen $n$',{'RMMSE($\varepsilon[t]$)' ''},[],[],{'ARMAX' 'ARX' 'BJ' 'OE'},'LegLoc','NorthEast');
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\Anwendungen\R51_NSFLTI_Cost.pdf','FigPosOffset',[0 0 10 0]);

% Estimate reduced Model
params.iopt.Hp = 1;	% Prediction horizon in #steps for identification
nA = 2;
nB = [1 3 1 1];
nC = 1;
M_nl = idModels.NsfPolyModel(nA,nB,ones(size(nB)),nC,params.mdl);
M_nl.factorize({'A' 'B'});
M_nl.identify(y,u2,params.iopt);
M_nl.updateSs('MinrealTol',1e-6); % Do PZ-Cancellation of StateSpace Object
E = cell2mat(M_nl.calcResiduals(y,u2,params.HpVal));
E = E(:); E = E(~isnan(E));
RMMSE = sqrt(mean(E.^2)); % RMMSE ERROR
T_1 = idModels.util.calcTimeconstant(M_nl.A.val(2:end),params.mdl.Ts)

% Plot BoxPlot and ACF
h1 = figure('color','w','Name','BoxAcf','Position',[100 100 800 300]); 
s(1)=subplot(1,2,1); col = lines(2);
M_arimax{1}.show('Boxplot',y,u1,'Hp',params.HpEst(end),'Axes',s(1),'Color',col(1,:)); legend('off');  
M_nl.show('Boxplot',y,u2,'Hp',params.HpEst(end),'Axes',s(1),'Color',col(2,:)); legend('off');ylim([-.35 .35])
s(2)=subplot(1,2,2); 
M_arimax{1}.show('Autocovariance',y,u1,'MaxLag',72,'Significance',.01,'Axes',s(2),'color',col(1,:)); legend('off');
M_nl.show('Autocovariance',y,u2,'MaxLag',72,'Significance',.01,'Axes',s(2),'color',col(2,:)); legend('off')
util.formatFigure(14);
%util.saveTightFigure(h1,'\\tsclient\D\Diss\Bilder\Anwendungen\R51_HammNl_BoxAcf.pdf','AxPosOffset',[-0.04 .06 0 0],'FigPosOffset',[0 0 -60 25])

% Plot Timeseries
col = lines(2);     % Colors for Plotting
ix_s = 6;           % Index of Dataset to be plotted 
fsz = 14;           % Fontsize
h2 = M_arimax{1}.show( 	'Prediction',y{ix_s},u1{ix_s},'Hp',params.HpVal,'TimeUnit','hours','FontSize',fsz,'ShowInputs',{1:2 5 6},...
                        'LeftAx',{[1 0] [1 1] []},'LineSpec',{{'-k' '-m'} {'-r'} {'-b'}},'HighlightOutliersByNsigma',[],'Color',col(1,:));
set(h2,'Position',[100 100 900 900])

figure(h2)
[~,yp] = M_nl.calcResiduals(y{ix_s},u2{ix_s},params.HpVal,Inf,true); 
subplot(h2.Children(8)); legend('off'); hold on; ylabel('$\vartheta_\mathrm{R}$  [$^{\circ}$C]','Interpreter','latex','FontSize',fsz) 
plot((0:length(yp)-1)*M_nl.Ts/60,squeeze(yp(:,1,end)),'Color',col(2,:),'DisplayName','Pr\"adiktion nichtlin. Modell mit Ausgangsr\"uckkopplung [$^{\circ}$C]');
ylim([23.9 25.6])

subplot(h2.Children(6)); legend('off'); ylabel('$\dot{Q}_\mathrm{Gs}$ [W/m$^2$]','Interpreter','latex','FontSize',fsz); yyaxis left; ylim([10 55]); 
hold on; plot((0:length(yp)-1)*M_nl.Ts/60,u2{ix_s}(:,8),'-r','DisplayName','$\vartheta_\mathrm{Vl,H}$ [$^{\circ}$C]'); 
hold on; plot((0:length(yp)-1)*M_nl.Ts/60,u2{ix_s}(:,10),'-b','DisplayName','$\vartheta_\mathrm{Vl,K}$ [$^{\circ}$C]'); hold on;
t1 = text(-2.5,3,'$\vartheta_\mathrm{Aul},$','Interpreter','latex','Color','k','FontSize',fsz,'Rotation',90);
t2 = text(t1.Position(1),t1.Position(2)+13.25,'$\vartheta_\mathrm{Vl,H},$','Interpreter','latex','Color','r','FontSize',fsz,'Rotation',90);
t3 = text(t2.Position(1),t2.Position(2)+15.5,'$\vartheta_\mathrm{Vl,K}$','Interpreter','latex','Color','b','FontSize',fsz,'Rotation',90);
t4 = text(t3.Position(1),t3.Position(2)+16.5,'[$^{\circ}$C]','Interpreter','latex','Color','k','FontSize',fsz,'Rotation',90);

s=subplot(h2.Children(4)); ylim([0 4]); legend('off');  ylabel('$\dot{Q}_\mathrm{H}$ [kW]','Interpreter','latex','FontSize',fsz); 
yyaxis right; hold on; ylabel('$\bar{H}_\mathrm{H}$','Interpreter','latex','FontSize',fsz); 
plot((0:length(yp)-1)*M_nl.Ts/60,mean(u2{ix_s}(:,5:7),2),'-k','DisplayName','$\bar{H}_\mathrm{H}$'); hold on;
set(s.YAxis(2),'TickLabelInterpreter','latex','Color','k');
set(s.YAxis(1),'Color','r');

s=subplot(h2.Children(2)); legend('off'); ylabel('$\dot{Q}_\mathrm{K}$ [kW]','Interpreter','latex','FontSize',fsz); 
yyaxis right; hold on; ylabel('$\bar{H}_\mathrm{K}$','Interpreter','latex','FontSize',fsz);  
plot((0:length(yp)-1)*M_nl.Ts/60,mean(u2{ix_s}(:,9),2),'-k','DisplayName','$\bar{H}_\mathrm{K}$'); hold on; ylim([0 1.1])
set(s.YAxis(2),'TickLabelInterpreter','latex','Color','k');
set(s.YAxis(1),'Color','b');
%util.saveTightFigure(h2,'\\tsclient\D\Diss\Bilder\Anwendungen\R51_Prediction.pdf','EqualX',1,'SubPlotYspace',.04,'FigPosOffset',[0 0 10 -210],'AxPosOffset',[0.01 -0.01 0 0]);

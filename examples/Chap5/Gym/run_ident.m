clc; clear; close all;
addpath(genpath('..\..\..'));

%% Load Data
params.Fname = 'Data_Gym.xls'; % Name of data file
[Info, Sheets] = xlsfinfo(params.Fname); % Load xlsinfo to get number of sheets
D = cellfun(@(sh) xlsread(params.Fname,sh),Sheets(2:end)','UniformOutput',false); % Load Data to tables
[~,ID] = xlsread(params.Fname,'Data1','B1:X1'); % Load column names        

%% Set Id Options
params.Ts = xlsread(params.Fname,'DataInfo','B3');      % Sampling Time of Model 
params.HpVal = 6*60/params.Ts;                         	% Prediction horizon in #steps for validation
params.Orders = 1:5;                                  	% Orders of Polynomials A,B,... for which identification will be performed
params.Type = 'armax';                                  % Parameterization for which identification will be performed
params.InitFromPreviosN = true;                         % Initialize Model of Order n with estimated Model with Order n-1
params.FactorizeG = true;                               % Identifiy model with factorized A(q),B(q),E(q)
params.fc = 1/(6*2*pi);                               	% Cutoff Frequency of highpass butterworth prefilter in 1/h. 
params.Seasonality = [7]*1440/params.Ts;            	% Seasonality of NoiseModel
params.SeasonalityOrder = 3;                            % Order of seasonal part
params.HpEst = [1 params.HpVal];                        % Prediction Horizon for estimation

% Options for Identification
params.iopt.TolX = 1e-8;                                % Stops if norm of change in parametervector is less then TolX
params.iopt.ForcePolesG = '';                           % Forces stable pos Real Poles to G 
params.iopt.MaxFunEval = 10e2;                        	% Max number of Cost function calls.  
params.iopt.MaxIter = Inf;                              % Max number of iterations during optimization.  
params.iopt.FunTolAbs = 1e-8;                           % Abort criterion of optimization. Abort if |F_i - F_i-1| <= FunTolAbs
params.iopt.StabilizationMethod = 'none';               % Stabilization of opt. predictor
params.iopt.Ic = 'backcast';                            % Initialization Method of residual filter 
params.iopt.IntegrateNoise = false;                   	% Force intergrator in noise model. 
params.iopt.InitMethod = 'arx';                         % Initialize Model with least squares
params.iopt.Solver = 'matlab_lsqnonlin';             	% Chosen optimization function. If no optimization toolbox installed you can use 'opti_levmar' from OPTI-Toolbox instead. (See https://inverseproblem.co.nz/OPTI/)
params.iopt.Algorithm = 'levenberg-marquardt';          % Chosen optimization algorithm (Remove this line if Solver = opti_levmar)
params.iopt.MultiStep = true;                           % Use MultiStep if Hp>1

%% HAMMERSTEIN MODEL: Extract Data and Identify
% Model Setup
params.col_in = [3 9]; % Outsidetemp., Heating power 
params.col_out = 1; % Roomtemperture 
params.mdl(1).Name = 'LTI';
params.mdl(1).Ts = params.Ts;
params.mdl(1).TimeUnit = 'minutes';
params.mdl(1).InputName = ID(params.col_in); 
params.mdl(1).OutputName  = ID(params.col_out);
params.mdl(1).InputNonlinearity(1).fun = '';        	% Identity
params.mdl(1).InputNonlinearity(1).parameters = [];
params.mdl(1).InputNonlinearity(1).free = []; 
params.mdl(1).InputNonlinearity(1).input_idx = 2;      % HeatingPow
params.mdl(1).InputNonlinearity(2).fun = '';           % Identity
params.mdl(1).InputNonlinearity(2).parameters = [];
params.mdl(1).InputNonlinearity(2).free = []; 
params.mdl(1).InputNonlinearity(2).input_idx = 1;      % OutsideTemp
params.mdl(1).PreFilter = cell(1,2); 
[params.mdl(1).PreFilter{:}] = butter(1,(params.fc/60)/((1/params.Ts)/2),'high');
params.mdl(1).PreFilter{1} = params.mdl(1).PreFilter{1}./params.mdl(1).PreFilter{1}(1);

% Extract Data
[y,u1] = cellfun(@(S) deal(S(1:end,params.col_out),S(1:end,params.col_in)),D,'UniformOutput',false); 

% % Do identification 
% for k = 1:length(params.HpEst)
%     params.iopt.Hp = params.HpEst(k);
%     [RES_lti{k},M_lti{k}] = idModels.test.sweepPar('Order',@idModels.NsfPolyModel.estModel,y,u1,'Order',params.Orders,...
%                             'Type',params.Type,'InitFromPreviousN',params.InitFromPreviosN,'MdlProp',params.mdl(1),...
%                            	'HpVal',params.HpVal,'FactorizeG',params.FactorizeG,'IdOpt',params.iopt,...
%                            	'SeasonalityC',params.Seasonality,'SeasonalityD',params.Seasonality,...
%                             'SeasonalityOrder',params.SeasonalityOrder);
% end

%% NONLINEAR MODEL (with supply temp. Measurement): Extract Data and Identify
% Model Setup
params.col_in = [3 10 8]; % Outsidetemp., Valve Pos, Supply Temp
params.col_out = 1; % Roomtemperture 
params.mdl(2).Name = 'NSF-LTI';
params.mdl(2).Ts =  params.Ts;
params.mdl(2).TimeUnit = 'minutes';
params.mdl(2).InputName = ID(params.col_in); 
params.mdl(2).OutputName  = ID(params.col_out);
params.mdl(2).InputNonlinearity(1).fun = 'idModels.func.fun_VdT'; 	% virt heating Power
params.mdl(2).InputNonlinearity(1).parameters = [1.1];           	% V*(Tsup-Troom)^p(1)
params.mdl(2).InputNonlinearity(1).free = [false]; 
params.mdl(2).InputNonlinearity(1).input_idx = [2 3];               % Valve Supply Temp   
params.mdl(2).InputNonlinearity(1).output_idx = 1;               	% Room Temp  
params.mdl(2).InputNonlinearity(2).fun = '';                        % Identity
params.mdl(2).InputNonlinearity(2).parameters = [];
params.mdl(2).InputNonlinearity(2).free = []; 
params.mdl(2).InputNonlinearity(2).input_idx = 1;                   % OutsideTemp
params.mdl(2).PreFilter = cell(1,2); 
[params.mdl(2).PreFilter{:}] = butter(1,(params.fc/60)/((1/params.Ts)/2),'high');
params.mdl(2).PreFilter{1} = params.mdl(2).PreFilter{1}./params.mdl(2).PreFilter{1}(1);

% Extract Data
[y,u2] = cellfun(@(S) deal(S(1:params.Ts/15:end,params.col_out),S(1:params.Ts/15:end,params.col_in)),D,'UniformOutput',false); 

% Do identification
% params.iopt.Hp = params.HpEst(1);
% [RES_nl,M_nl] = idModels.test.sweepPar({'Order'},@idModels.NsfPolyModel.estModel,y,u2,'Order',params.Orders,...
%                         'Type',params.Type,'InitFromPreviousN',params.InitFromPreviosN,'MdlProp',params.mdl(2),...
%                         'HpVal',params.HpVal,'FactorizeG',params.FactorizeG,'IdOpt',params.iopt,...
%                         'SeasonalityC',params.Seasonality,'SeasonalityD',params.Seasonality,...
%                         'SeasonalityOrder',params.SeasonalityOrder);

%% Estimate specific Model and compare model with and without seasonality
% Params
nA = 1;
nB = [1 1];
nC = 0;
nD = 0;
nS = 1;  

% Identify Seasonal Model
ixD = false(1,max([nD+1 params.Seasonality+1+(nS-1)])); ixD(2:nD+1) = true;
for ns = 1:length(params.Seasonality)
    ixD(params.Seasonality(ns) + 1 + (-nS+1:nS-1)) = true;
end
ixC = false(1,max([nC+1 params.Seasonality+1+(nS-1)])); ixC(2:nC+1) = true;
for ns = 1:length(params.Seasonality)
    ixC(params.Seasonality(ns) + 1 + (-nS+1:nS-1)) = true;
end
M1 = idModels.NsfPolyModel(nA,nB,ones(size(nB)),length(ixC)-1,length(ixD)-1,params.mdl(1));
M1.C.val([false ~ixC(2:end)]) = 0;
M1.C.free = ixC;
M1.D.val([false ~ixD(2:end)]) = 0;
M1.D.free = ixD;
M1.factorize({'A' 'B'});
M1.identify(y,u1,params.iopt);
idModels.util.calcTimeconstant(M1.A.val(2),M1.Ts)/60

M2 = idModels.NsfPolyModel(nA,nB,ones(size(nB)),length(ixC)-1,length(ixD)-1,params.mdl(2));
M2.C.val([false ~ixC(2:end)]) = 0;
M2.C.free = ixC;
M2.D.val([false ~ixD(2:end)]) = 0;
M2.D.free = ixD;
M2.factorize({'A' 'B'});
M2.identify(y,u2,params.iopt);

%% Plot BoxPlot and ACF
h1 = figure('color','w','Name','BoxAcf','Position',[100 100 800 300]); 
s(1)=subplot(1,2,1); col = lines(2);
M1.show('Boxplot',y,u1,'Hp',params.HpVal,'Axes',s(1),'Color',col(1,:),'DeltaK',2); legend('off');  
M2.show('Boxplot',y,u2,'Hp',params.HpVal,'Axes',s(1),'Color',col(2,:),'DeltaK',2,'LineStyle','--'); legend('off'); %ylim([-.65 .65]);
s(1).XTick = 0:4:12; s(1).XTickLabel = arrayfun(@num2str,0:2:6,'UniformOutput',false); xlabel('Pr\"adiktionshorizont in h');
s(2)=subplot(1,2,2); 
M1.show('Autocovariance',y,u1,'MaxLag',4*24*8,'Significance',.01,'Axes',s(2),'color',col(1,:)); legend('off');
M2.show('Autocovariance',y,u2,'MaxLag',4*24*8,'Significance',.01,'Axes',s(2),'color',col(2,:),'LineStyle','--'); legend('off'); ylim([-.25 1.0]*10e-4)
s(2).XTick = [0 2 4 6 8]*96; s(2).XTickLabel = arrayfun(@num2str,[0 2 4 6 8],'UniformOutput',false); xlabel('Lag $\tau$ in d');
util.formatFigure(14);
%util.saveTightFigure(h1,'\\tsclient\D\Diss\Bilder\Anwendungen\Gym_LtiNl_BoxAcf.pdf','AxPosOffset',[-0.04 .07 0 0],'FigPosOffset',[0 0 -60 35])

%% Prediction
ix_set = 1;      	% Index of Dataset to be plotted   
col = lines(2);     % line colors
fsz = 12;           % Fontsize
h2 = M1.show(   'Prediction',y{ix_set},u1{ix_set},'Hp',params.HpVal,'TimeUnit','days','ShowInputs',{1 2},...
             	'LeftAx',{1 1},'LineSpec',{{'-k'} {'-r'}},'Color',col(1,:),'FontSize',fsz); 
set(h2,'Position',[100 100 900 900])
[~,yp] = M2.calcResiduals(y{ix_set},u2{ix_set},params.HpVal,Inf,true); 

% s1
subplot(h2.Children(6)); legend('off'); hold on; ylabel('$\vartheta_\mathrm{R}$  [$^{\circ}$C]','Interpreter','latex','FontSize',fsz) 
plot((0:length(yp)-1)*M2.Ts/1440,squeeze(yp(:,1,end)),'Color',col(2,:),'LineStyle','-');
l = legend({'Messung' 'LTI Modell' 'Nichtlin. Modell'},'Interpreter','latex','FontSize',fsz,'Location','southwest');
%set(l,'Position',[0.11,0.71,l.Position(3),0.095])
ylim([17.0 22]);

% s2
s=subplot(h2.Children(4)); legend('off'); hold on; ylabel('$\vartheta_\mathrm{Aul}$  [$^{\circ}$C]','Interpreter','latex','FontSize',fsz);
yyaxis right; stairs((0:length(yp)-1)*M2.Ts/1440,u2{ix_set}(:,3),'-r'); ylabel('$\vartheta_\mathrm{Vl,H}$  [$^{\circ}$C]','Interpreter','latex','FontSize',fsz);
set(s.YAxis(2),'TickLabelInterpreter','latex','Color','r');

% s3
s=subplot(h2.Children(2)); legend('off'); hold on; ylabel('$\dot{Q}_\mathrm{H}$ [kW]','Interpreter','latex','FontSize',fsz);
yyaxis right; stairs((0:length(yp)-1)*M2.Ts/1440,u2{ix_set}(:,2),'-k'); ylabel('$H_P$','Interpreter','latex','FontSize',fsz);
set(s.YAxis(2),'TickLabelInterpreter','latex','Color','k'); xlim([10 40]);
s.XTickLabel = arrayfun(@num2str,0:5:30,'UniformOutput',false);
%util.saveTightFigure(h2,'\\tsclient\D\Diss\Bilder\Anwendungen\Gym_Prediction.pdf','EqualX',1,'SubPlotYspace',.07,'FigPosOffset',[0 0 0 -230],'AxPosOffset',[0.0 0.01 0 0]);
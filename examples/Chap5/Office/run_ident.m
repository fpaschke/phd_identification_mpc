clc; clear; close all;
addpath('..\..\..')

%% Load Data
params.Fname = 'Data_Office.xls'; % Name of data file
[Info, Sheets] = xlsfinfo(params.Fname); % Load xlsinfo to get number of sheets
D = cellfun(@(sh) xlsread(params.Fname,sh),Sheets(2:end)','UniformOutput',false); % Load Data to tables
[~,ID] = xlsread(params.Fname,'Data1','B1:X1'); % Load column names        

%% Set Id Options  
params.Ts = xlsread(params.Fname,'DataInfo','B3');      % Sampling Time of Model 
params.Win = xlsread(params.Fname,'DataInfo','B6:B7');  % Orientation of Windows (wrt North and to a hor. plane)
params.HpVal = 3*60/params.Ts;                         	% Prediction horizon in #steps for validation
params.Orders = 1:2;                                 	% Orders of Polynomials A,B,... for which identification will be performed
params.Type = 'armax';                                	% Parameterization for which identification will be performed
params.InitFromPreviosN = false;                      	% Initialize Model of Order n with estimated Model with Order n-1
params.FactorizeG = true;                               % Identifiy model with factorized A(q),B(q),E(q)
params.fc = 1/(6*2*pi);                                 % Cutoff Frequency of highpass butterworth prefilter in 1/h. 
params.Seasonality = {[] [7]*1440/params.Ts};        	% Seasonality of NoiseModel
params.SeasonalityOrder = 2;                            % Order of Seasonal Part

% Options for Identification
params.iopt.TolX = 1e-8;                                % Stops if norm of change in parametervector is less then TolX
params.iopt.ForcePolesG = '';                       	% Forces stable pos Real Poles to G if 'sp'
params.iopt.MaxFunEval = 10e2;                        	% Max number of Cost function calls.  
params.iopt.MaxIter = Inf;                              % Max number of iterations during optimization.  
params.iopt.FunTolAbs = 1e-8;                           % Abort criterion of optimization. Abort if |F_i - F_i-1| <= FunTolAbs
params.iopt.Ic = 'backcast';                            % Initialization Method of residual filter 
params.iopt.IntegrateNoise = false;                   	% Force intergrator in noise model. 
params.iopt.InitMethod = 'arx';                         % Initialize Model with least squares
params.iopt.Solver = 'matlab_lsqnonlin';             	% Chosen optimization function. If no optimization toolbox installed you can use 'opti_levmar' from OPTI-Toolbox instead. (See https://inverseproblem.co.nz/OPTI/)
params.iopt.Algorithm = 'trust-region-reflective';  	% Chosen optimization algorithm (Remove this line if Solver = opti_levmar)

%% NONLINEAR MODEL: Extract Data and Identify
% Model Setup
params.col_in = [2 3 4 5 6 7]; % Outsidetemp., global radiation on horizontal plane, sun height, sun azimuth, heating valve pos., heating supply temp
params.col_out = 1; % Roomtemperture 
params.mdl = struct();
params.mdl.Name = 'NSF-LTI';
params.mdl.Ts =  params.Ts;
params.mdl.TimeUnit = 'minutes';
params.mdl.InputName = ID(params.col_in); 
params.mdl.OutputName  = ID(params.col_out);
params.mdl.InputNonlinearity(1).fun = 'idModels.func.fun_VdT'; % "virt. Heatingpower": Vdot * (Tsupply_heat - Troom)^p1 
params.mdl.InputNonlinearity(1).parameters = [1.1]; % p1 
params.mdl.InputNonlinearity(1).free = [false]; 
params.mdl.InputNonlinearity(1).input_idx = 5:6; % Valve, SupplyTemp
params.mdl.InputNonlinearity(1).output_idx = 1; % Feedback of room temp
params.mdl.InputNonlinearity(2).fun = ''; % Identity function for outside temp.
params.mdl.InputNonlinearity(2).parameters = [];
params.mdl.InputNonlinearity(2).free = []; 
params.mdl.InputNonlinearity(2).input_idx = 1; % OutsideTemp
params.mdl.InputNonlinearity(3).fun = 'idModels.func.fun_rad'; % calculates radiation on inclined plane and multiply by logistic function to consider shading by trees and buildings in front of the window
params.mdl.InputNonlinearity(3).parameters = [params.Win' 2]; % orientation of windows wrt N in degrees, orientation of windows wrt to horizontal plane in degrees, estimated shading angle, growth rate/steepness of logistic function
params.mdl.InputNonlinearity(3).free = [false false false]; 
params.mdl.InputNonlinearity(3).input_idx = [2 4 3]; % Radiation, Azimuth, Height    
params.mdl.PreFilter = cell(1,2); 
[params.mdl.PreFilter{:}] = butter(1,(params.fc/60)/((1/params.Ts)/2),'high');
params.mdl.PreFilter{1} = params.mdl.PreFilter{1}./params.mdl.PreFilter{1}(1);

% Extract Data
[y,u] = cellfun(@(S) deal(S(1:params.Ts/15:end,params.col_out),S(1:params.Ts/15:end,params.col_in)),D,'UniformOutput',false); 

%% Do identification and plot RMMSE
for k = 1:length(params.Seasonality)
    [RES_nl{k},M_nl{k}] = idModels.test.sweepPar({'Order'},@idModels.NsfPolyModel.estModel,y,u,'Order',params.Orders,...
                            'Type',params.Type,'InitFromPreviousN',params.InitFromPreviosN,'MdlProp',params.mdl,...
                           	'HpVal',params.HpVal,'FactorizeG',params.FactorizeG,'IdOpt',params.iopt,...
                           	'SeasonalityC',params.Seasonality{k},'SeasonalityD',params.Seasonality{k},...
                            'SeasonalityOrder',params.SeasonalityOrder);
end

% Plot RMMSE
for k = 1:length(params.Seasonality) 
  	[mat(:,k), lab] = util.tab2mat(RES_nl{k},1); 
end
figure('color','w','Name','NSF-LTI','Position',[200 200 500 280]); bar(mat); grid on; ylim([.2 .28])
util.formatFigure(14,'$n$',{'RMMSE($\varepsilon[t]$)' ''},[],[],{{'ARMAX' 'S-ARARMAX'}},'LegLoc','NorthEast');
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\Anwendungen\B7_Cost.pdf','FigPosOffset',[0 0 5 -5],'AxPosOffset',[-0.00 -0.02 0 0.0])

%% Estimate specific Model and compare model with and without seasonality
% Params
nA = 2;
nB = [1 1 1]; 
nC = 1;
nD = 0;
nS = 2;

% Identify Seasonal Model
ixD = false(1,max([nD+1 params.Seasonality{2}+1+(nS-1)])); ixD(2:nD+1) = true;
for ns = 1:length(params.Seasonality{2})
    ixD(params.Seasonality{2}(ns) + 1 + (-nS+1:nS-1)) = true;
end
ixC = false(1,max([nC+1 params.Seasonality{2}+1+(nS-1)])); ixC(2:nC+1) = true;
for ns = 1:length(params.Seasonality{2})
    ixC(params.Seasonality{2}(ns) + 1 + (-nS+1:nS-1)) = true;
end
Ms = idModels.NsfPolyModel(nA,nB,ones(size(nB)),length(ixC)-1,length(ixD)-1,params.mdl);
ixC(end) = false;
Ms.C.val([false ~ixC(2:end)]) = 0;
Ms.C.free = ixC;
ixD(end-1:end) = false;
Ms.D.val([false ~ixD(2:end)]) = 0;
Ms.D.free = ixD;
Ms.factorize({'A' 'B'});
Ms.identify(y,u,params.iopt);
E = cell2mat(Ms.calcResiduals(y,u,params.HpVal));
E = E(:); E = E(~isnan(E));
RMMSEs = sqrt(mean(E.^2)); % RMMSE ERROR
T = idModels.util.calcTimeconstant(Ms.A.val(2:end),params.mdl.Ts)/60
Ms.printParameters;

% IdentifyModel  with no seasonality
Mns = idModels.NsfPolyModel(nA,nB,ones(size(nB)),nC,params.mdl);
Mns.factorize({'A' 'B'});
Mns.identify(y,u,params.iopt);
E = cell2mat(Mns.calcResiduals(y,u,params.HpVal));
E = E(:); E = E(~isnan(E));
RMMSEns = sqrt(mean(E.^2)); % RMMSE ERROR

% Plot BoxPlot and ACF
h1 = figure('color','w','Name','BoxAcf','Position',[100 100 900 350]); 
s(1)=subplot(1,2,1); col = lines(2);
Mns.show('Boxplot',y,u,'Hp',params.HpVal,'Axes',s(1),'Color',col(1,:)); legend('off');  
Ms.show('Boxplot',y,u,'Hp',params.HpVal,'Axes',s(1),'Color',col(2,:)); legend('off'); ylim([-.65 .65])
s(1).XAxis.TickValues = [0 4 8 12]; s(1).XAxis.TickLabels = {'0' '1' '2' '3'}; s(1).XLabel.String = 'Pr\"adiktionshorizont $k$ [h]';
s(2)=subplot(1,2,2); 
Mns.show('Autocovariance',y,u,'MaxLag',700,'Significance',.01,'Axes',s(2),'color',col(1,:)); legend('off');
Ms.show('Autocovariance',y,u,'MaxLag',700,'Significance',.01,'Axes',s(2),'color',col(2,:)); legend('off'); ylim([-3 100]*1e-4);
s(2).XAxis.TickValues = (0:7)*96; s(2).XAxis.TickLabels = {'0' '1' '2' '3' '4' '5' '6' '7'}; s(2).XLabel.String = 'Lag $\tau$ [d]';
util.formatFigure(14);
%util.saveTightFigure(h1,'\\tsclient\D\Diss\Bilder\Anwendungen\B7_BoxAcf.pdf','AxPosOffset',[-0.04 .04 0 0],'FigPosOffset',[0 0 -60 15])

% Plot Bode Diagramm 
Ms.showBode('Legend',{},'ShowFrequencyWeighting',false,'ShowH',false);
set(gcf,'Position',[100 100 1000 300])
util.formatFigure(14); 
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\Anwendungen\B7_Bode.pdf','FigPosOffset',[0 0 -60 30],'AxPosOffset',[-.07 0.07 0 0],'SubPlotXSpace',-.01);
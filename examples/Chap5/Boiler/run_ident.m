clc; clear; close all;
addpath('..\..\..')

%% Load Data
params.Fname = 'Data_Gasboiler.xls';        % Name of data file
[~, Sheets] = xlsfinfo(params.Fname);    % Load xlsinfo to get number of sheets
D = cellfun(@(sh) xlsread(params.Fname,sh),Sheets(2:end)','UniformOutput',false); % Load Data to tables
[~,ID] = xlsread(params.Fname,'Data1','B1:E1'); % Load column names        

%% Set Id Options  
params.Ts = xlsread(params.Fname,'DataInfo','B3');      % Sampling Time of Model 
params.Pnom = xlsread(params.Fname,'DataInfo','B4');    % Nominal Heating Power [kW]
params.HpVal = 5/params.Ts;                         	% Prediction horizon in #steps for validation
params.UseConstraints = false;                        	% Use Constraints if true

% Options for Identification
params.iopt.TolX = 1e-8;                                % Stops if norm of change in parametervector is less then TolX
params.iopt.ForcePolesG = '';                           % Forces stable pos Real Poles to G if 'sp'
params.iopt.MaxFunEval = 5e2;                        	% Max number of Cost function calls.  
params.iopt.MaxIter = 5e2;                              % Max number of iterations during optimization.  
params.iopt.FunTolAbs = 1e-8;                           % Abort criterion of optimization. Abort if |F_i - F_i-1| <= FunTolAbs
params.iopt.StabilizationMethod = 'none';               % Stabilization of opt. predictor
params.iopt.Ic = 'ls';                                  % Initialization Method of residual filter 
params.iopt.IntegrateNoise = false;                   	% Force intergrator in noise model. 
params.iopt.EstimateOutputOffset = false;             	% Estimates OutputOffset if true
params.iopt.InitMethod = 'arx';                         % Initialize Model with least squares
params.iopt.Solver = 'matlab_lsqnonlin';             	% Chosen optimization function. If no optimization toolbox installed you can use 'opti_levmar' from OPTI-Toolbox instead. (See https://inverseproblem.co.nz/OPTI/)
params.iopt.Algorithm = 'levenberg-marquardt';          % Chosen optimization algorithm (Remove this line if Solver = opti_levmar)

%% MODEL OF BOILER
% Model Setup
params.col_in = [4 2 3];    % Mod, Tret, VolFlow
params.col_out = 1;         % Tsup
params.mdl = struct();
params.mdl.Name = 'NSF-LTI';
params.mdl.Ts =  params.Ts;
params.mdl.TimeUnit = 'minutes';
params.mdl.InputName = ID(params.col_in); 
params.mdl.OutputName  = ID(params.col_out);

% Setup InputNonlinearity
params.mdl.InputNonlinearity(1).fun = 'idModels.func.fun_VdT'; % V*(Tret-Tsup): Convection Into Boiler
params.mdl.InputNonlinearity(1).parameters = 1;
params.mdl.InputNonlinearity(1).free = 0;
params.mdl.InputNonlinearity(1).input_idx = [3 2];
params.mdl.InputNonlinearity(1).output_idx = 1;
params.mdl.InputNonlinearity(2).fun = 'idModels.func.fun_boiler_eff'; % Heat Gain due combustion
params.mdl.InputNonlinearity(2).output_idx = [];
params.mdl.InputNonlinearity(2).parameters = [params.Pnom 100 50 .1 .11]';
params.mdl.InputNonlinearity(2).free = [0 0 1 1 0];
params.mdl.InputNonlinearity(2).input_idx = [1 2];
if params.UseConstraints 
    params.mdl.InputNonlinearity(2).A = [0 0 -1  0  0   %p(3)>=30
                                         0 0  1  0  0   %p(3)<=70
                                         0 0  0 -1  0   %p(4)>=1e-3
                                         0 0  0  1  0   %p(4)<=1                                        
                                         0 0  0  0 -1   %p(5)>=0
                                         0 0  0  0  1]; %p(5)<=0.2
    params.mdl.InputNonlinearity(2).b = [-30 70  -1e-3 1 0 .2]';
end
    
% Extract Data
[y,u] = cellfun(@(S) deal(S(1:end,params.col_out),S(1:end,params.col_in)),D,'UniformOutput',false); 

M = idModels.NsfPolyModel(2,2,1,2,0,0,'Ts',params.mdl.Ts,'TimeUnit',params.mdl.TimeUnit,...
    'OutputName',params.mdl.OutputName,'InputName',params.mdl.InputName,'InputNonlinearity',params.mdl.InputNonlinearity);
M.factorize({'A'}); % Factorize A (and B): Not necessary for 1st order Models
M.A.val(2:end) = [1 .9]; M.A.free(2:end) = [0 1]; 
M.C.val = conv([1 -.95],[1 -.95]); M.C.free(2:end) = false;
M.identify(y,u,params.iopt);
[~,yp,misc] = M.calcResiduals(y,u,params.HpVal,Inf,true); 
misc.Rmmse 
M.printParameters
M.showBode;
idModels.util.calcTimeconstant(M.A.val(end),1)

%% ABB XXX
Sets = [68 120 169]; Ns = length(Sets);
N = cellfun(@length,yp)-1;

figure('color','w','Position',[200 200 900 500]); 
for ns = 1:length(Sets)
    s(1,ns) = subplot(2,3,ns); stairs(0:N(Sets(ns)),y{Sets(ns)},'-r'); hold on; stairs(0:N(Sets(ns)),u{Sets(ns)}(:,2),'-b'); hold on; plot(0:N(Sets(ns)),yp{Sets(ns)}(:,1,end),'-k'); grid on; if ns == 1; ylabel('Temperatur [$^{\circ}$C]'); end% legend({'Vorlauftemp.' 'R\"ucklauftemp.' 'Pr\"ad. Vorlauftemp.'},'Location','southeast');
    s(2,ns) = subplot(2,3,ns+Ns); stairs(0:N(Sets(ns)),u{Sets(ns)}(:,1),'-k'); xlabel('Zeit $t$ [min]'); if ns == 1; ylabel('Modulation [\%]'); end; ylim ([0 100]);  hold on; yyaxis right; stairs(0:N(Sets(ns)),u{Sets(ns)}(:,3)/3.6,'-m'); ylim ([0 12]/3.6); if ns == 3; ylabel('Volumenstrom [l/s]','Color','m'); end; grid on;  s(2,ns).YAxis(2).Color = 'm';
    set(gca,'Position',get(gca,'Position')+[0 0.05 0 0]);
end
util.formatFigure(14);
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\Anwendungen\ThermePred.pdf','XShift',-.05,'FigPosOffset',[0 0 -70 -45],'AxPosOffset',[0 -.05 0 0])

%% ABB XXX
f = M.showBode('ShowH',true,'ShowPhase',1,'Labels','Sparse'); 
set(f,'Position',[200 200 1000 280]); util.formatFigure(14);
subplot(1,3,1); set(gca,'XTick',[.01 1]);
subplot(1,3,2); set(gca,'XTick',[.01 1]);
subplot(1,3,3); set(gca,'XTick',[.01 1]);
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\Anwendungen\ThermeBode.pdf','XShift',-.075,'FigPosOffset',[0 0 -70 30],'AxPosOffset',[.0 .075 0 0],'SubplotXSpace',-.01)

%% ABB XXX
M.InputMin(2) = 25; M.InputMax(2) = 70;
f = M.showInputNonlinearity(2,[100 NaN],'Confidence',.99,'Color','k'); xlim([25 70]); ylim([250 290]); legend off;
set(f,'Position',[200 200 400 280]); 
util.formatFigure(14,'R\"ucklauftempertur $\vartheta_\mathrm{Rl}$ [$^{\circ}$C]','Heizleistung $\dot{Q}_\mathrm{K}$ [kW]');
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\Anwendungen\ThermeQb.pdf','FigPosOffset',[0 0 0 5]);

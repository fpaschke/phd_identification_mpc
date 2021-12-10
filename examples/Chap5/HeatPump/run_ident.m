clc; clear; close all;
addpath('..\..\..')

%% Load Data
params.Fname = 'Data_HeatPump.xls';   	% Name of data file
[~, Sheets] = xlsfinfo(params.Fname); 	% Load xlsinfo to get number of sheets
D = cellfun(@(sh) xlsread(params.Fname,sh),Sheets(2:end)','UniformOutput',false); % Load Data to tables
[~,ID] = xlsread(params.Fname,'Data1','B1:I1'); % Load column names        

%% Set Id Options
% Model Parameters
params.Ts = xlsread(params.Fname,'DataInfo','B3');      % Sampling Time of Model 
params.P0 = 50;                                         % Nominal Heating Power of Heatpump if one compressor is active T0 [kW]
params.T0 = [10 35];                                    % Nominal Temperature PrimSup/SecSup [C]
params.Cop0 = 5.6;                                      % Nominal COP at Tnom
params.P0Comp = 0;                                      % Compressor constant loss [kW]
params.EtaComp = .9;                                    % Compressor efficiency [0...1] 0.9 means that 90% of El. Power is converted to Heat that will be available for heating
params.HpVal = 5/params.Ts;                         	% Prediction horizon in #steps for validation
params.UseConstraints = true;                        	% Use Constraints for parameters of InputNonlinearities
params.Order = 1;                                       % Polynomial Orders 

% Options for Identification
params.iopt.TolX = 1e-12;                               % Stops if norm of change in parametervector is less then TolX
params.iopt.ForcePolesG = 'sp';                         % Forces stable pos Real Poles to G 
params.iopt.MaxFunEval = 1e3;                        	% Max number of Cost function calls.  
params.iopt.MaxIter = 1e3;                              % Max number of iterations during optimization.  
params.iopt.FunTolAbs = 1e-8;                           % Abort criterion of optimization. Abort if |F_i - F_i-1| <= FunTolAbs
params.iopt.StabilizationMethod = 'none';               % Stabilization of opt. predictor
params.iopt.Ic = 'ls';                                  % Initialization Method of residual filter 
params.iopt.IntegrateNoise = false(2,1);                % Force intergrator in noise model. 
params.iopt.EstimateOutputOffset = false(2,1);         	% Estimates OutputOffset if true
params.iopt.InitMethod = 'arx';                         % Initialize Model with least squares
params.iopt.Solver = 'matlab_fmincon';              	% Chosen optimization function. If no optimization toolbox installed you can use 'opti_levmar' from OPTI-Toolbox instead. (See https://inverseproblem.co.nz/OPTI/)
params.iopt.Algorithm = 'sqp';                          % Chosen optimization algorithm (Remove this line if OPTI Toolbox is used)
params.iopt.Hp = 1;                                     % Prediction Horizon
params.iopt.MultiStep = false;                        	% Use Multistep Criterion if 1

%% MODEL SETUP 
params.col_in = [4 1 6 5];      % TsecRet,TprimSup,Vsec,Vprim
params.col_out = [2 3];         % TsupSec,TretPrim 
params.mdl = struct();
params.mdl.Name = 'HeatpumpModel';
params.mdl.Ts =  params.Ts;
params.mdl.TimeUnit = 'minutes';
params.mdl.InputName = ['Ansteuerung' ID(params.col_in)]; %Input 1: Compressor1 + Compressor2 (see data extraction)
params.mdl.OutputName  = ID(params.col_out);

% SECONDARY InputNonlinearities
% Convection Into Condenser
params.mdl.InputNonlinearity(1).fun = 'idModels.func.fun_VdT';  %Vsec*(TsecRet - TsupSec)
params.mdl.InputNonlinearity(1).parameters = 1;
params.mdl.InputNonlinearity(1).free = 0;
params.mdl.InputNonlinearity(1).input_idx = [4 2]; %Vsec, TsecRet 
params.mdl.InputNonlinearity(1).output_idx = 1; %TsecSup

% Heat Gain due condensation
params.mdl.InputNonlinearity(2).fun = 'idModels.func.fun_HPpow_sec'; %x/p(1)*(p(2) + p(5)*dTss + p(6)*dTps + p(7)*dTss*dTps + p(8)*dTss^2 + p(9)*dTps^2) with dTss = Tss - p(3) and dTps = Tps - p(4)
params.mdl.InputNonlinearity(2).parameters = [1 params.P0 params.T0(2) params.T0(1) -0.003*params.P0 0.03*params.P0 0 0 0];
if params.UseConstraints
    params.mdl.InputNonlinearity(2).A = [	0 0 0 0 1  0 0 0  0;    % p(5)<0
                                        	0 0 0 0 0 -1 0 0  0;    % p(6)>0
                                            0 0 0 0 0  0 0 0 -1];   % p(9)>0 
    params.mdl.InputNonlinearity(2).b = [0 0 0]';
end
params.mdl.InputNonlinearity(2).free = [0 0 0 0 1 1 0 0 0]; %x0 P0 Tsv0 Tpv0 dTsv dTpv dTsv*dTpv dTsv^2 dTpv^2
params.mdl.InputNonlinearity(2).input_idx = [1 3]; %X, TprimSup  
params.mdl.InputNonlinearity(2).output_idx = 1; %TsupSec 

% PRIMARY InputNonlinearities
% Convection Into Evaporator
params.mdl.InputNonlinearity(3).fun = 'idModels.func.fun_VdT'; %Vprim*(TprimSup-TprimRet)
params.mdl.InputNonlinearity(3).parameters = 1;
params.mdl.InputNonlinearity(3).free = 0;
params.mdl.InputNonlinearity(3).input_idx = [5 3]; %Vprim, TprimSup 
params.mdl.InputNonlinearity(3).output_idx = 2; 

% Heat Loss due evaporation
dT0 = params.T0(2)-params.T0(1);
params.mdl.InputNonlinearity(4).fun = 'idModels.func.fun_HPpow_prim'; %Psec*(COP-etaComp)/COP - P0Comp 
params.mdl.InputNonlinearity(4).parameters = [params.P0Comp params.EtaComp params.Cop0 dT0 -0.15 7.34e-4 ...
                   arrayfun(@(np) ['2_' num2str(np)],1:9,'UniformOutput',false)];
params.mdl.InputNonlinearity(4).free = [0 0 0 0 1 1];
if params.UseConstraints
    Nc = 3;
    params.mdl.InputNonlinearity(4).A = [zeros(Nc,4) ones(Nc,1) 2*linspace(20-dT0,60-dT0,Nc)'];
    params.mdl.InputNonlinearity(4).b = zeros(Nc,1);
end
params.mdl.InputNonlinearity(4).input_idx = [1 3]; %X,TprimSup  
params.mdl.InputNonlinearity(4).output_idx = 1; %TsupSec 

%% Extract Data and identify
[y,u] = cellfun(@(S) deal(S(:,params.col_out),[sum(S(:,7:8),2) S(:,params.col_in)]),D,'UniformOutput',false); 

params.Order = 1;
Mprop = util.struct2namevaluepairs(params.mdl);
M = idModels.NsfPolyModel(params.Order*eye(2),params.Order*[1 1 0 0; 0 0 1 1],ones(2,4),params.Order*ones(2,1),Mprop{:});
M.factorize({'A' 'B'});
M.A(1,1).val(2:end) = [1]; M.A(1,1).free(2:end) = [0]; 
M.A(2,2).val(2:end) = [1]; M.A(2,2).free(2:end) = [0];
M.C(1).val = conv([1 -.95],[1 ]); M.C(1).free(2:end) = 0; 
M.C(2).val = conv([1 -.95],[1 ]); M.C(2).free(2:end) = 0;
M.identify(y,u,params.iopt);

M.show('Prediction',y(1:10),u(1:10),'Hp',params.HpVal,'ObserverInitSamples',10,'ShowInputs',{2 3 1 [4 5]},'LineSpec',{{'-k' '--r'  '-b'} {'-k' '--b'  '-y'} {'-k'} {'--r' '--b'}})

%% ABB XXX
figure('Color','w','Position',[200 200 1000 320]); N = 100; alpha = .99;
dT = linspace(20,50,N)';
Tprim = 6:3:15; Nt = length(Tprim);
for k = 1:Nt
    s1 = subplot(1,3,1); M.showInputNonlinearity(4,[1 Tprim(k)],NaN,'Axes',s1,'color','k','Confidence',alpha);
    text(33+2*k,33+4*k,[num2str(Tprim(k)) '$^{\circ}$C'],'Interpreter','latex','FontSize',14);
    s2 = subplot(1,3,2); M.showInputNonlinearity(2,[1 Tprim(k)],NaN,'Axes',s2,'color','k','Confidence',alpha);
    text(32+2*k,37+4.7*k,[num2str(Tprim(k)) '$^{\circ}$C'],'Interpreter','latex','FontSize',14);
end
[f,df] = idModels.func.fun_poly_normalized(dT,[M.InputNonlinearity(4).parameters{3:6}]);
dF = zeros(N,length(M.Info.Popt));
dF(:,end-sum(M.InputNonlinearity(4).free(3:6))+1:end) = df(:,M.InputNonlinearity(4).free(3:6));
stdF = sqrt(diag(dF*M.Info.CovP*dF'));
confF = util.norminv((alpha+1)/2,zeros(N,1),stdF);      
s3 = subplot(1,3,3); plot(dT,f,'Color','k'); hold on; plot(dT,f+[confF -confF],'LineStyle','--','Color','k'); grid on;
util.formatFigure(14,{'$\vartheta_\mathrm{sek,Vl}$ [C$^{\circ}$]' '$\vartheta_\mathrm{sek,Vl}$ [C$^{\circ}$]' '$\vartheta_\mathrm{sek,Vl}-\vartheta_\mathrm{prim,Vl}$ [K]'},...
    {'Verdampferleistung $\dot{Q}_\mathrm{prim}$ [kW]' 'Verfl\"ussigerleistung $\dot{Q}_\mathrm{sek}$ [kW]' 'Leistungszahl $\varepsilon$'});
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\Anwendungen\Wpf.pdf','XShift',-.075,'FigPosOffset',[0 0 -90 30],'AxPosOffset',[.0 .075 0 0],'SubplotXSpace',-.01)

%% ABB XXX
figure('Color','w','Position',[200 200 1000 450]);
for ny = 1:M.OutputDimension
    s = subplot(2,3,1+3*(ny-1)); M.showBode('Axes',s,'Idx',[ny 1+2*(ny-1)]); if ny == 1; xlabel(''); end; set(gca,'XTick',[.01 1]); yyaxis left; ylabel('Amplitude $A(\omega)$ [dB]'); yyaxis right; ylabel(''); 
    s = subplot(2,3,2+3*(ny-1)); M.showBode('Axes',s,'Idx',[ny 2+2*(ny-1)]); if ny == 1; xlabel(''); end; set(gca,'XTick',[.01 1]); yyaxis right; ylabel('');
    s = subplot(2,3,3+3*(ny-1)); M.showBode('Axes',s,'Idx',[ny]); if ny == 1; xlabel(''); end; set(gca,'XTick',[.01 1]);
end
util.formatFigure(14);
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\Anwendungen\WpBode.pdf','XShift',-.07,'FigPosOffset',[0 0 -95 0],'AxPosOffset',[0 .01 0 0],'SubplotYSpace',-.00,'SubplotXSpace',-.0)

%% ABB XXX
Sets = [1 2 19]; Ns = length(Sets);
Y = y(Sets); U = u(Sets);
[~,yp] = M.calcResiduals(Y,U,'Hp',params.HpVal);
N = cellfun(@length,yp)-1;
figure('color','w','Position',[200 200 900 500]); 
for ns = 1:Ns
    s(1,ns) = subplot(2,3,ns); 
    stairs(0:N(ns),Y{ns}(:,1),'-r'); hold on; 
    stairs(0:N(ns),U{ns}(:,2),'-b'); hold on; 
    plot(0:N(ns),yp{ns}(:,1,end),'-k'); grid on; 
    xlim([0 20]);
    if ns == 1; ylabel('Temperatur [$^{\circ}$C]'); end 
    yyaxis right; l = stairs(0:N(ns),U{ns}(:,4)/3.6,'-m'); ylim([0 4]);
    s(1,ns).YAxis(2).Color = l.Color;
    if ns == length(N); ylabel('$\dot{V}_\mathrm{sek}$ [l/s]'); end 
	s(2,ns) = subplot(2,3,ns+Ns); 
    stairs(0:N(ns),U{ns}(:,3),'-r'); hold on; 
    stairs(0:N(ns),Y{ns}(:,2),'-b'); hold on; 
    plot(0:N(ns),yp{ns}(:,2,end),'-k'); grid on; 
    xlim([0 20]);
    if ns == 1; ylabel('Temperatur [$^{\circ}$C]'); end 
    yyaxis right; l = stairs(0:N(ns),U{ns}(:,5)/3.6,'-m'); ylim([0 7]);
    s(2,ns).YAxis(2).Color = l.Color;
    if ns == length(N); ylabel('$\dot{V}_\mathrm{prim}$ [l/s]'); end
    xlabel('$t$ in min');
end
util.formatFigure(14);
%util.saveTightFigure(gcf,'D:\Diss\Bilder\Anwendungen\WpPred.pdf','XShift',-.05,'FigPosOffset',[0 0 -70 -20],'AxPosOffset',[0 .0 0 0])

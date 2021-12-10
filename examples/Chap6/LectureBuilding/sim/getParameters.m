function P = getParameters(varargin)
% Returns default parameter struct
%
% P = getParameters(varargin)
% P [struct]:	Simulation, Controller and Identification Parameters
% varargin: Optional name value pairs
% GENERAL PARAMETERS:
% T_s [integer]: Sampling Time in min.
% N_w [integer]: Number of simulation [weeks]
% IDENTIFICATION PARAMETERS
% Hp_val [integer]: Time Horizon for model validation in h.


p = inputParser(); 
addParameter(p,'T_s',5,@(x) x>=1 && mod(x,1)==0);
addParameter(p,'Hp_val',48); 
addParameter(p,'T_s_ident',60); 
addParameter(p,'T_s_mpc',60); 
addParameter(p,'N_w',52); 
parse(p,varargin{:})

%% Simulation and Simulink Model Parameters
rng(0);                                     % For Reproducability of N_pers
P.Tstart = '2021-1-1';                      % Assumed StartDate of simulation (year is irrelavant)
P.DateFormat = 'dd.mm HH:MM';               % Format for Plotting timeseries
P.N_w = p.Results.N_w;                      % Number of simulation [weeks]
P.T_s = p.Results.T_s;                      % Sampling time in [min]
P.T_sens = 1;                           	% Timeconstant of Temperature Sensors [min]
P.Q_int = .1;                           	% Quantization Interval of Temperature Sensors [K]
P.Eta_fan = 0.1;                          	% Power Coefficient of AHU-Fan [kW/(kg/s)^2]. Pel = 2*eta*m^2 where m^2 
P.P_fan = @(m,eta) deal(2* eta*m.^2,4*eta*m);% Computes el. power consuption [kW] of fans (2 fans) and derivative dP/dm 
P.Cost_heat = .1;                         	% Cost of Heating/Cooling Power [€/kW]        
P.Cost_el = .3;                             % Cost of Electrical Power [€/kW]        
P.Eta_heatrecovery = .75;               	% Average Efficiency of Heat Recovery
P.C_air = 1.005;                        	% Heat Capacity of Air [kJ/kgK]     
P.GeoPos = [52.47 13.40];                   % Geoposition of building (latitude and longitude) [degrees]
P.Solver = 'ode45';                         % Solver for FMU Integration
P.RelTol = '1e-6';                          % Relative tolerance for integration
P.MaxStep = num2str(P.T_s*60);              % Max timestep of integration [s]
P.WeatherFile = 'DEU_Berlin.103840_IWEC.mos'; % File for WeatherData to be read
P.CtrlFun = 'ctrl_pi';                      % Function name for controller to be used (see: .\subfun\ctrl)

%% Input Signal Parameters
P.Sig.N_pers.Schedule = [7.0 8.5; 9.0 10.5; 11.0 12.5;% Occupancy Schedule for Mo-Fr [hours of day]
                         13.5 15.0; 15.5 17.0; 17.5 19]+1;            
P.Sig.N_pers.Mean = round([300 350 250 200]'.* ...  	% Mean of Number of people in each zone [#Number of people]
                    rand(4,5,size(P.Sig.N_pers.Schedule,1))); %dim1: RoomNumber, dim2: day (Mo-Fr), dim3: hour  
P.Sig.N_pers.T_r = 15;                         	% RiseTime [min]
P.Sig.N_pers.Std = 0.05;                      	% Relative Standard dev. of number of people (relative to Mean) [1]
P.Sig.T_ref.DayTime = [P.Sig.N_pers.Schedule(1,1) ... 
                       P.Sig.N_pers.Schedule(end,2)];	% Start of day- and night-time during weekdays [hours of day]
P.Sig.T_ref.T_set_lower = [15 21].*ones(4,1); 	% min Zone Setpoints (nighttime/daytime) [C]
P.Sig.T_ref.T_set_upper = [30 23].*ones(4,1); 	% max Zone Setpoints (nighttime/daytime) [C]
% Upper bound in summer will be computed based on maximum outside temp. (lin. interp)
P.Sig.T_ref.T_set_summer = 26;                  % Upper bound summer temp.
P.Sig.T_ref.T_out_summer = [26 32];          	% Outsidetemp summer     
P.Sig.T_ref.PrbsAmplitude = 0.*ones(4,1);       % Add PRBS with specified Amplitude to setpoint [K] (Can be used to generate more informative data for identification)
P.Sig.T_ref.PrbsInterval = 1.*ones(4,1);        % New random value each x [h]
P.Sig.T_ref.PrbsPhaseShift = 0;                 % Switching point of signal i will be x [hours] after switching point of signal i-1

% Signals for Identification
r = p.Results.T_s_ident/60;
P.Sig.T_sup.Range = [18 26];                 	% Range of random Signal [C]
P.Sig.T_sup.StepSize = 1;                       % StepSize of random Signal [K]
P.Sig.T_sup.SwitchInterval = 48*r;          	% New random value each x [h]     
P.Sig.T_sup.TimeSpan = [0 P.N_w*7];           	% Start and End of Signal [day]
P.Sig.m_sup = repmat(struct(...
                'Range',[0 6],...               % Range of random Signal [m3/s]
                'StepSize',1,...              	% StepSize of random Signal [m3/s]
                'SwitchInterval',24*r,...     	% New random value each x [h]       
                'TimeSpan',[0 P.N_w*7]),1,4); 	% Start and End of Signal [day]
P.Sig.m_sup(2).TimeSpan = [0 P.N_w*7]+6/24*r; 
P.Sig.m_sup(3).TimeSpan = [0 P.N_w*7]+12/24*r;
P.Sig.m_sup(4).TimeSpan = [0 P.N_w*7]+18/24*r;
P.Sig.m_sup(3).Range = [0 4];
P.Sig.m_sup(4).Range = [0 4];

%% General Controller Parameters (Necessary for PI and MPC Control)
P.Ctrl.T_sup_day = [18 26];                     % Supply Temps. of AHU during day (cooling/heating-mode)
P.Ctrl.T_sup_night = [14 30];                   % Supply Temps. of AHU during night (cooling/heating-mode)
P.Ctrl.m_max = [6 6 4 4]';                      % Maximum Massflowrate to each zone [kg/s]
P.Ctrl.T_day = [P.Sig.N_pers.Schedule(1,1)     	% Start of day- and night-time during weekdays 
             	P.Sig.N_pers.Schedule(end,2)];	% (detrmines max and min supply temperatures) [hours of day]              

%% PI Controller parameters
P.Ctrl.PI.T_pre = 0.5*[1 1 1 1];            	% Preheating Time in [h]
P.Ctrl.PI.T_s = P.T_s;                          % Sampling time (used to compute integral output) [min]
P.Ctrl.PI.Hyst = [P.Sig.T_ref.T_set_lower(:,2) ....
                  P.Sig.T_ref.T_set_upper(:,2)];% Switching points for avarage hall temperature for cooling and heating mode [C]
P.Ctrl.PI.K_p = 1*[1 1 1 1]';                   % Propotional gains of PI-Controlers [(kg/s)/K)]
P.Ctrl.PI.K_i = 6*P.Ctrl.PI.T_s/60*[1 1 1 1]';  % Integral Gains of PI-Controlles [(kg/s^2)/(K)]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ]

if false
    M = load(['.\resource\Models_Ts' num2str(P.Ctrl.PI.T_s) '.mat']);
    M.RoomMod_SO.defactorize();
    for nz = 1:length(P.Ctrl.PI.T_pre)
        Gi = tf(M.RoomMod_SO.B(nz,1+2*nz).val,M.RoomMod_SO.A(nz,nz).val,M.RoomMod_SO.Ts,'Variable','z^-1');
        dT = 5; % Linearization Point (Tsup-Troom)
        [P.Ctrl.PI.K_p, P.Ctrl.PI.K_i] = getPI(dT*Gi,'DiscMth','zoh',...
                            'Ts',P.Ctrl.PI.T_s*60,'Method','tsum');
    end
end
%% MPC Controller parameters
P.Ctrl.MPC.T_s = p.Results.T_s_mpc;                                                                
P.Ctrl.MPC.ModelFile = ['.\resource\Models_Ts' num2str(P.Ctrl.MPC.T_s) '.mat'];   % Path to Model File
P.Ctrl.MPC.MpcActivationTime = P.Tstart;            % Prior MpcActivationTime the PI Controller will be used (Useful for debugging). 
P.Ctrl.MPC.UseDeltaU = false;                       % Uses differeneces of (Unsupported yet!)
P.Ctrl.MPC.UseSeasonalModel = false;                % If true then a seasonal model for occupancy modelling will be used (currently not available if P.Ctrl.MPC.UseMimoModel = true)
P.Ctrl.MPC.UseMimoModel = false;                   	% If true MIMO Model will be used. If false MISO Models will be used. See: run_ident.m
P.Ctrl.MPC.UseWeatherPredictionModel = false;   	% If true weather predictions will be used. If false exact weather data will be used.
P.Ctrl.MPC.UseSparse = false;                       % If true sparse state tranbsition matrices will be used.
P.Ctrl.MPC.Wconst = 1e3;                            % Constraint violation factor for cost function V = ... + P.Ctrl.MPC.Wconst*Hy*eps^2 (eps shifts constraint bounds) 
if p.Results.T_s_mpc == 120
    P.Ctrl.MPC.Hy = 1:12*7;
    P.Ctrl.MPC.Hu = P.Ctrl.MPC.Hy-1;
else
%     P.Ctrl.MPC.Hu = [0:23 24:2:(4*24-2)];        	% Control horizon during night and weekend. Should start with 0 (= current timestep)
%     P.Ctrl.MPC.Hy = [1:24 26:2:4*24];            	% Prediction horizon during night and weekend.    
  	P.Ctrl.MPC.Hu = 0:95;                        	% Control horizon during night and weekend. Should start with 0 (= current timestep)
    P.Ctrl.MPC.Hy = 1:96;                           % Prediction horizon during night and weekend.    
end

% Solver Options: see fmincon options                 
P.Ctrl.MPC.SolOpt.Algorithm = 'sqp';                % Alg to be used (Active-set, sqp or interior-point)
P.Ctrl.MPC.SolOpt.SpecifyObjectiveGradient = true; 	% Supply Gradient for CostFun
P.Ctrl.MPC.SolOpt.SpecifyConstraintGradient = true; % Supply Gradient for Constraints
P.Ctrl.MPC.SolOpt.CheckGradients = false;         	% Check analytical Gradient
P.Ctrl.MPC.SolOpt.Display = 'off';                	% Display iterations
P.Ctrl.MPC.SolOpt.ConstraintTolerance = 1e-2;   	% Termination Tolerance for satisfying constraints
P.Ctrl.MPC.SolOpt.FiniteDifferenceType = 'central'; % Only used if CheckGradients=true
P.Ctrl.MPC.SolOpt.OptimalityTolerance = 1e-6;       % Termination Tolerance on first order optimality
P.Ctrl.MPC.SolOpt.MaxFunctionEvaluations = 5e2;     % Max Number of function calls

%% Identification Patrameters
P.Ident.MIMO = false;                       % Identify MIMO ARMAX Model if true
P.Ident.T_s = p.Results.T_s_ident;       	% Sampling time of identified model [min]
P.Ident.Orders = 1:3;                    	% Orders of Polynomials A,B,... for which identification will be performed
P.Ident.Type = 'armax';                 	% Model Structure that will be used for identification     
P.Ident.FactorizeG = true;               	% Identifiy model with factorized A(q),B(q),E(q)
P.Ident.MinRealTol = [];                    % Perform model order reduction of ss object using minreal(ss,tol)
P.Ident.fc = [];                            % Cutoff Frequency of highpass butterworth prefilter in 1/h. (1/(3*2*pi) corresponds to Timconstant of 3 hours);                  
P.Ident.HpVal = p.Results.Hp_val*60/P.Ident.T_s;	% Prediction horizon in #steps for validation
P.Ident.InitFromPrevious = false;        	% Use previous estimated Model for Parameter initialization
P.Ident.WinDir = [270 90 270 90];           % Orientation of windows wrt to North (W,E,W,E)

% Options for Identification 
P.Ident.Iopt.TolX = 1e-10;                	% Stops if norm of change in parametervector is less then TolX
P.Ident.Iopt.ForcePolesG = '';          	% Forces stable pos Real Poles to G if 'sp'
P.Ident.Iopt.MaxIter = 2e3;                 % Max number of Iterations 
P.Ident.Iopt.MaxFunEval = 2e3;          	% Max number of Cost function calls.  
P.Ident.Iopt.FunTolAbs = 1e-10;          	% Abort criterion of optimization. Abort if |F_i - F_i-1| <= FunTolAbs
P.Ident.Iopt.Ic = 'backcast';             	% Initialization Method of residual filter 
P.Ident.Iopt.IntegrateNoise = false;     	% Force intergrator in noise model. 
P.Ident.Iopt.InitMethod = 'arx';         	% Initialize Model with least squares
P.Ident.Iopt.Solver = 'matlab_lsqnonlin'; 	% Chosen optimization function. If no optimization toolbox installed you can use 'opti_levmar' from OPTI-Toolbox instead. (See https://inverseproblem.co.nz/OPTI/)
P.Ident.Iopt.Algorithm = 'levenberg-marquardt';    	% Chosen optimization algorithm (Remove this line if Solver = opti_levmar)

%% Checks
% for nz = 1:length(P.Ctrl.MPC.Models.RoomMod)
%     if P.T_s ~= P.Ctrl.MPC.Models.RoomMod{nz}.Ts/60
%         P.Ctrl.MPC.Models.RoomMod{nz}.resample(P.T_s*60,'zoh')
%     end
%     %assert(P.T_s==P.Ctrl.MPC.Models.RoomMod{nz}.Ts/60,'Models do not have specified sampling time! Consider resampling first!');
% end
end
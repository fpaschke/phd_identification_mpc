function P = getParameters(varargin)
% Returns default parameter struct
%
% P = getParameters(varargin)
% P [struct]:	Simulation, Controller and Identification Parameters
% varargin: Optional name value pairs
% GENERAL PARAMETERS:
% T_s [integer]: Sampling Time in min.
% PeriodicOccupancy [logical]: generates periodic occupancy if true Def. false
% IDENTIFICATION PARAMETERS
% Hp_val [integer]: Time Horizon for model validation in h.

p = inputParser(); 
addParameter(p,'T_s',5,@(x) x>=1 && mod(x,1)==0);
addParameter(p,'Hp_val',6); 
addParameter(p,'PeriodicOccupancy',false,@(x) x==1 || x==0); 
parse(p,varargin{:})

%% Simulation and Simulink Model Parameters
rng(0);                                         % For Reproducability of N_pers
P.Tstart = '2021-1-1';                          % Assumed StartDate of simulation (year is irrelavant)
P.DateFormat = 'dd.mm HH:MM';                   % Format for Plotting timeseries
P.N_w = 52;                                     % Number of simulation [weeks]
P.T_s = p.Results.T_s;                          % Sampling time in [min]
P.T_sens = 1;                                   % Timeconstant of Temperature Sensor [min]
P.Q_int = 1e-1;                                 % Quantization Interval of Temperature Sensors [K]
P.GeoPos = [50.827 	12.921];                	% Geoposition of room (latitude and longitude) [degrees]
P.Solver = 'ode45';                           	% Solver for FMU Integration
P.RelTol = '1e-6';                              % Relative tolerance for integration of FMU (to low values can cause errors)
P.MaxStep = num2str(P.T_s*60);                  % Max timestep of integration [s]
P.WeatherFile = 'DEU_Berlin.103840_IWEC.mos'; 	% File for WeatherData to be read
P.CtrlFun = 'ctrl_pi';                          % Function name for controller to be used (see: .\subfun\ctrl)
P.CoolDays = [120 273];                         % Days of cooling (120=1.5.XXXX 273=1.10.XXXX)
if ~p.Results.PeriodicOccupancy
    P.ModelFile = ['.\resource\Armax_Ts' num2str(P.T_s) '_Hp15m.mat']; % Path to model file that will be used for MPC and PI Control design
else
	P.ModelFile = ['.\resource\S-Armax_Ts' num2str(P.T_s) '_Hp15m.mat']; % Path to model file that will be used for MPC and PI Control design
end
%% Input Signal Parameters
P.Sig.N_pers.Periodic = p.Results.PeriodicOccupancy; % If true then a periodic occupancy pattern will be generated, otherwise an random pattern           
if ~P.Sig.N_pers.Periodic
    P.Sig.N_pers.NumPeop = 8;                 	% Max/Mean Number of people [#Number of people] (random/periodic)
    P.Sig.N_pers.Schedule = 9:2:18;            	% Occupancy Schedule for Mo-Fr [hours of day]' (random)  
else
    P.Sig.N_pers.NumPeop = 6;                 	% Max/Mean Number of people [#Number of people] (random/periodic)
    P.Sig.N_pers.Schedule = [9 13; 14 18];          % Occupancy Schedule for Mo-Fr [hours of day]' (periodic)            
end
P.Sig.N_pers.StdPeop = .2;                    	% relative Std. Deviation of people [% Npeop]
P.Sig.N_pers.RiseTime = 30;                     % RiseTime [min]

P.Sig.T_ref.DayTime = [ P.Sig.N_pers.Schedule(1)-2 ... 
                        P.Sig.N_pers.Schedule(end)+2];                    % Start of day- and night-time during weekdays [hours of day]
P.Sig.T_ref.T_set_lower = [16 21];           	% min Zone Setpoints (nighttime/daytime) [C]
P.Sig.T_ref.T_set_upper = [30 23];              % max Zone Setpoints (nighttime/daytime) [C]
% Upper bound in summer will be computed based on maximum outside temp. (lin. interp)
P.Sig.T_ref.T_set_summer = 26;                  % Upper bound summer temp.
P.Sig.T_ref.T_out_summer = [26 32];          	% Outsidetemp summer     

% Signals for Identification
P.Sig.H_flo.Range = [0 1];                      % Range of random Signal [1]
P.Sig.H_flo.StepSize = .1;                      % StepSize of random Signal [1]
P.Sig.H_flo.SwitchInterval = 12;                % New random value each x [h] 
P.Sig.H_flo.TimeSpan = [0 7*12];              	% New random value each x [h] 
P.Sig.H_win.Range = [0 1];                      % Range of random Signal [1]
P.Sig.H_win.StepSize = .25;                    	% StepSize of random Signal [1]
P.Sig.H_win.SwitchInterval = 48;                % New random value each x [h] 
P.Sig.H_win.TimeSpan = [0 7*12];             	% New random value each x [h] 

%% Controller parameters for "standard" control  
% Floor heating
P.Ctrl.Std.Flo.T_pre = 6;                    	% Preheating Time in [h]
P.Ctrl.Std.Flo.T_s = P.T_s;                   	% Sampling time of floor heating controller [min]
P.Ctrl.Std.Flo.H_max = 1;                     	% Maximum Valve Position [1]
% Controller parameters that have been chosen empirically 
P.Ctrl.Std.Flo.K_p = 1;                        	% Propotional gain of PI-Controler [1/K)]
P.Ctrl.Std.Flo.K_i = .5/60*P.Ctrl.Std.Flo.T_s; 	% Integral gain of PI-Controller [1/(K)]

% Design PI Controller using Empirical Method (Tsummen Regel, Betragsopt, Chien/Hrones/Reswick)
if true
    M = load(P.ModelFile);Sys = tf(M.RoomModFlo.Ss); n = M.RoomModFlo.InputNonlinearity(3).parameters{2}; m = M.RoomModFlo.InputNonlinearity(3).parameters{1};
    H = .5; % Linearization point Valve
    dT = 15; % Linearization Point (Tsup-Troom)
    [P.Ctrl.Std.Flo.K_p, P.Ctrl.Std.Flo.K_i] = getPI(m*(H)^(m-1)*dT^n*Sys(3),'DiscMth','zoh',...
                        'Ts',P.Ctrl.Std.Flo.T_s*60,'Method','absopt');
end

% Shading
P.Ctrl.Std.Win.Mode = 'bb';                     % Determines weather shading is used as control input
                                                % 'off': blinds will be open all the time (no control input)
                                                % 'time': blinds will closed in coolmode (P.CoolDays) and opened otherwise  
                                                % 'time_unocc': same as time but if occupied then blinds are opened
                                                % 'bb': uses bang bang control to determine shading position
P.Ctrl.Std.Win.T_s = 180;                      	% Sampling Time [min]. If 60min the shading will be changed only at full hour                              
P.Ctrl.Std.Win.BbHyst = [P.Sig.T_ref.T_set_lower(2) P.Sig.T_ref.T_set_upper(2)]; %
                                                

%% MPC Controller parameters 
% Sampling time is the sampling time of simulation and needs to match model sampling time (model specified by P.ModelFile)
P.Ctrl.MPC.MpcActivationTime = P.Tstart;      	% Prior MpcActivationTime the PI Controller will be used (Useful for debugging). 
P.Ctrl.MPC.UseWeatherPredictionModel = false; 	% If true weather predictions will be used. If false exact weather data will be used.
P.Ctrl.MPC.Wy_day = 0;                        	% Tracking error weight day
P.Ctrl.MPC.Wy_noday = 0;                      	% Tracking error weight no day
P.Ctrl.MPC.UseOutputConstraints = 'quad';      	% Choose the following:
                                                % 'off' deactivates Output Constraints
                                                % 'lin'/'quad' enables output constraints and uses linear/quadratic penalization method ie V = ... + Wconst* gamma (gamma^2) where gamma shifts comfort bounds
P.Ctrl.MPC.Wconst = 1000;                    	% Constraint violation factor if UseOutputConstraints = true for cost function V = ... + P.Ctrl.MPC.Wconst*eps (eps shifts upper and lower constraint bounds) 
P.Ctrl.MPC.Hy_noday = 12*4;                 	% Prediction steps during night and weekend.
P.Ctrl.MPC.Hy_day = 12*4;                    	% Prediction steps during day.
P.Ctrl.MPC.T_day = P.Sig.T_ref.DayTime;       	% Start of day- and night-time during weekdays (which prediction/control horizon will be used) [hours of day]            

% Floor heating
P.Ctrl.MPC.Flo.Wu_day = 1e-1;                	% Weight for dQv (virt. heating power) during day
P.Ctrl.MPC.Flo.Wu_noday = 1e-1;              	% Weight for dQv (virt. heating power) during night and weekend
P.Ctrl.MPC.Flo.Hu_noday = P.Ctrl.MPC.Hy_noday;  % Control hor during night and weekend.
P.Ctrl.MPC.Flo.Hu_day = P.Ctrl.MPC.Hy_day;  	% Control hor during daytime.
% NOTE: It is good to choose Hu = Hy since the computation of maxmimum vitual power is approximated for k>Hu.


% Shading
P.Ctrl.MPC.Win.Mode = 'on';                     % Determines shading mode:
                                                % 'on': blind pos will be treated as control input. 
                                                % 'on_unocc' blind pos will be treated as control input if room is unoccupied. 
                                                % 'closed'/'opened' leaves blinds in specified position
P.Ctrl.MPC.Win.T_s = P.Ctrl.Std.Win.T_s;      	% Sampling Time [min]. If 60min the shading will be changed only at full hour
P.Ctrl.MPC.Win.Wu_day = 1e-3;               	% Weight for dHv (shading)
P.Ctrl.MPC.Win.Wu_noday = 1e-3;               	% Weight for dHv (shading)
P.Ctrl.MPC.Win.Hu_noday = 4;                 	% Control hor during night and weekend. (Wrt P.Ctrl.MPC.Win.T_s)
P.Ctrl.MPC.Win.Hu_day = 4;                    	% Control hor during daytime. (Wrt P.Ctrl.MPC.Win.T_s)

% Solver
P.Ctrl.MPC.SolOpt.Display = 'off';             	% Display iterations
P.Ctrl.MPC.SolOpt.Algorithm = 'active-set';    	% quadprog algorithm. 

%% Identification Patrameters
P.Ident.Orders = 1:5;                   	% Orders of Polynomials A,B,... for which identification will be performed
P.Ident.Type = 'armax';                   	% Model Structure that will be used for identification     
P.Ident.FactorizeG = true;               	% Identifiy model with factorized A(q),B(q),E(q)
P.Ident.fc = [];                            % Cutoff Frequency of highpass butterworth prefilter in 1/h. (1/(3*2*pi) corresponds to Timconstant of 3 hours);                  
P.Ident.HpVal = p.Results.Hp_val*60/P.T_s;	% Prediction horizon in #steps for validation
P.Ident.InitFromPrevious = true;        	% Use previous estimated Model for Parameter initialization
P.Ident.WinDir = 270;                       % Orientation of window wrt to North (W=270)
if P.Sig.N_pers.Periodic
    P.Ident.Seasonality = 24*7*60/P.T_s;   
else
    P.Ident.Seasonality = [];    
end
P.Ident.SeasonalityOrder = 'ByOrder';		% Order of seasonal part (pos int.), if 'ByOrder' then it will be increased with P.Ident.Orders 

% Options for Identification 
P.Ident.Iopt.TolX = 1e-8;                   % Stops if norm of change in parametervector is less then TolX
P.Ident.Iopt.ForcePolesG = 'sp';          	% Forces stable pos Real Poles to G if 'sp'
P.Ident.Iopt.ForceZerosG = '';              % Forces pos Real Zeros to G if 'p'
P.Ident.Iopt.MaxFunEval = 5e2;          	% Max number of Cost function calls.  
P.Ident.Iopt.MaxIterations = 5e2;        	% Max number of Iterations  
P.Ident.Iopt.FunTolAbs = 1e-6;          	% Abort criterion of optimization. Abort if |F_i - F_i-1| <= FunTolAbs
P.Ident.Iopt.Ic = 'backcast';              	% Initialization Method of residual filter 
P.Ident.Iopt.IntegrateNoise = false;     	% Force intergrator in noise model. 
P.Ident.Iopt.InitMethod = 'arx';         	% Initialize Model with least squares
P.Ident.Iopt.Solver = 'matlab_lsqnonlin'; 	% Chosen optimization function. If no optimization toolbox installed you can use 'opti_levmar' from OPTI-Toolbox instead. (See https://inverseproblem.co.nz/OPTI/)
P.Ident.Iopt.Algorithm = 'trust-region-reflective';    	% Chosen optimization algorithm (Remove this line if Solver = opti_levmar)

end
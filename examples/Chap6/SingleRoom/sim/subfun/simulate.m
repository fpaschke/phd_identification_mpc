function [SimOut, t_sim] = simulate(sys,P)
%Simulates Simulink Model sys and calculates energy and comfort.   
%The necessary parameters and Inputs are assumed to be loaded into
%workspace. Execute getSimIn and getSimPar prior calling this function.

% Clear all local function workspaces
clear ctrl_fun ctrl_mpc ctrl_pi;

% Gen Input Signals and write them to base workspace
In = getSimIn(P);                       
assignin('base','In',In)

% Load Simulink Model with FMU
load_system(sys);                         

% Set Simulation Parameters and create Bus for weather Data
fmudialog.createBusType([sys '/SingleRoomFMU/FMU']);	
set_param(sys,'Solver',P.Solver,'RelTol',P.RelTol,'MaxStep',P.MaxStep,...
        'StartTime',num2str(In.Time(1)),'StopTime',num2str(In.Time(1)+P.N_w*7*60*24*60),'SimulationCommand','update');
                         
% Simulate
ts = tic; 
SimOut = sim(sys); 
t_sim = toc(ts)

% If Open Loop Simulation
if ~isprop(SimOut,'H_flo')  
    SimOut.H_flo = In.H_flo;
    SimOut.H_win = In.H_win;
end

% Calculate Energies and Comfort
Ph = SimOut.P_flo.Data; Ph(Ph<0) = 0;
Pc = SimOut.P_flo.Data; Pc(Pc>0) = 0;
SimOut.E_heat_flo = timeseries(cumsum(P.T_s*Ph/60),SimOut.T_air.Time,'Name','Heating energy'); % [kWh]         
SimOut.E_cool_flo = timeseries(-cumsum(P.T_s*Pc/60),SimOut.T_air.Time,'Name','Cooling energy'); % [kWh]         

N = length(SimOut.T_air.Data);
e_lower = In.T_ref_lower.Data(1:N) - SimOut.T_air.Data;    	% Lower Control Error [K]        	
e_lower = e_lower.*(e_lower>0);            	% Only lower violations 
e_upper = SimOut.T_air.Data - In.T_ref_upper.Data(1:N);   	% Upper Control Error [K]        	
e_upper = e_upper.*(e_upper>0);            	% Only upper violations 

%ix_cooling = SimOut.T_air.Time>=P.CoolDays(1)*86400 & SimOut.T_air.Time<P.CoolDays(2)*86400;
SimOut.Comfort_lin = timeseries(cumsum(e_lower + e_upper)*P.T_s/60,SimOut.T_air.Time,'Name','Linear Comfort Measure'); % Linear Comfort measure [Kh]  
SimOut.In = In;                             % Add simulation inputs to SimOut
SimOut.Parameters = P;                      % Add parameters to SimOut            

% Set Timeproperties 
SimOut = setTimeProp(SimOut,[P.Tstart(1:4) '-1-1'],P.DateFormat);
SimOut.T_air.DataInfo.Interpolation = 'linear';
end


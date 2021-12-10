function [SimOut, t_sim] = simulate(sys,P)
%Simulates Simulink Model sys and calculates energy and comfort.   
%The necessary parameters and Inputs are assumed to be loaded into
%workspace. Execute getSimIn and getSimPar prior calling this function.

% Data to
global Uopt;                    % Saves optimal control input trajectories          
global Yopt;                    % Saves optimal control output trajectories
global Xobs;                    % Saves Observed state

% Clear all local function workspaces
clear ctrl_fun ctrl_mpc ctrl_pi;

% Gen Input Signals and write them to base workspace
In = getSimIn(P);                       
assignin('base','In',In)

% Load Simulink Model with FMU
load_system(sys);                         

% Set Simulation Parameters and create Bus for weather Data
fmudialog.createBusType([sys '/LectureBuildingFMU/FMU']);	
set_param(sys,'Solver',P.Solver,'RelTol',P.RelTol,'MaxStep',P.MaxStep,...
        'StartTime',num2str(In.Time(1)),'StopTime',num2str(In.Time(1)+P.N_w*7*60*24*60),'SimulationCommand','update');
                         
% Simulate
ts = tic; 
SimOut = sim(sys); 
t_sim = toc(ts)

% If Open Loop Simulation
if ~isprop(SimOut,'m_sup')  
    SimOut.m_sup = In.m_sup;
    SimOut.T_sup = In.T_sup;
end

% Calculate Energies and Comfort
SimOut.M_ahu = timeseries(sum(SimOut.m_sup.Data,2),SimOut.T_air.Time,'Name','AHU Massflowrate'); 
N = length(SimOut.M_ahu.Data); Nz = size(In.T_ref_lower.Data,2);
[P_el,~] = P.P_fan(SimOut.M_ahu.Data,P.Eta_fan);
SimOut.P_el = timeseries(P_el,SimOut.T_air.Time,'Name','AHU el. power'); %[kW]
SimOut.E_el = timeseries(cumsum(P.T_s*P_el)/60,SimOut.T_air.Time,'Name','AHU el. energy consump'); % [kWh]
SimOut.P_heat = timeseries((1-P.Eta_heatrecovery)*P.C_air*SimOut.M_ahu.Data.*(SimOut.T_sup.Data - ...
                    interp1(In.T_out.Time,In.T_out.Data,SimOut.T_sup.Time)),SimOut.T_air.Time,'Name','AHU heating/cooling power'); % [kW]
SimOut.E_heat = timeseries(cumsum(P.T_s*SimOut.P_heat.Data.*(SimOut.P_heat.Data>0))/60,SimOut.T_air.Time,'Name','AHU heating energy'); % [kWh]         
SimOut.E_cool = timeseries(-cumsum(P.T_s*SimOut.P_heat.Data.*(SimOut.P_heat.Data<0))/60,SimOut.T_air.Time,'Name','AHU cooling energy'); % [kWh]           
SimOut.Cost_el = P.Cost_el*SimOut.E_el;                  	% Cost for el. En. [€]
SimOut.Cost_heat = P.Cost_heat*SimOut.E_heat;           	% Cost for heat. En. [€]
SimOut.Cost_cool = P.Cost_heat*SimOut.E_cool;           	% Cost for cool. En. [€]
SimOut.Cost = SimOut.Cost_el + SimOut.Cost_heat + SimOut.Cost_cool; % Total Cost [€]
SimOut.e_lower = In.T_ref_lower.Data(1:N,:) - SimOut.T_air.Data(:,1:Nz);    % Lower Control Error during Occ [K]        	
SimOut.e_lower = SimOut.e_lower.*(SimOut.e_lower>0);                            % Only lower violations 
SimOut.e_lower_nz = sum(SimOut.e_lower,1);
SimOut.e_upper = SimOut.T_air.Data(:,1:Nz) - In.T_ref_upper.Data(1:N,:);    % Upper Control Error during Occ [K]      	
SimOut.e_upper = SimOut.e_upper.*(SimOut.e_upper>0);                            % Only upper violations 
SimOut.e_upper_nz = sum(SimOut.e_upper,1);
SimOut.Comfort_lin = timeseries(cumsum(sum(SimOut.e_lower + SimOut.e_upper,2))*P.T_s/60,SimOut.T_air.Time,'Name','LinearComfortMeasure');    	% Linear Comfort measure [Kh]  
SimOut.In = In;                                                                 % Add simulation inputs to SimOut
SimOut.Parameters = P;
SimOut.Uopt = Uopt;
SimOut.Yopt = Yopt;
SimOut.Uobs = Xobs;

% Set Timeproperties 
SimOut = setTimeProp(SimOut,[P.Tstart(1:4) '-1-1'],P.DateFormat);
SimOut.T_air.DataInfo.Interpolation = 'linear';
SimOut.Parameters = P;
end


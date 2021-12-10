function [T_sup, m_sup, Cool_mode, t_sol] = ctrl_pi(t,T_air)
%Implements a MIMO PI-Controller for maintaining T_air within the bounds T_ref.
%Each of the PI Controllers computes m_sup and can be in heating or cooling state. 
%The state determines if lower or upper setpoint is used (T_ref(1) or T_ref(2)).  
%The state switches using bang-bang-scheme (CtrlPar.Hyst).
%The supply temperature T_sup is also determined using bang-bang scheme. If
%at least 3 PI controllers are in cooling state then the AHU switches into
%cooling state as well which means that CtrlPar.T_sup(1,1) or CtrlPar.T_sup(2,1) will be
%used as the setpoint (determined by time). If at least 3 PI controllers
%are in heating state then AHU controller switches back into heating mode.
%
% [T_sup, m_sup, Cool_mode, t_sol] = ctrl_pi(t,T_air)
% T_sup [scalar double]:    Supply Temp in [C]
% m_sup [4 x 1 double]:     Massflow rate in [kg/s]
% Cool_mode [4 x 1 double]:	State of PI Controllers [0/1]
% t [scalar double]:        time in [s] (used to determine weather day-/nighttime or weekend)
% T_air [4 x 1 double]:     Measured Room Temperatures [C]

tic;
persistent Par;        	% Parameter struct 
persistent x;                   % State variables of Integral part of PI controllers
persistent state_coolmode;      % State variables of Bang Bang Controller of rooms
persistent state_coolmode_AHU; 	% State of AHU (cooling/heating)
persistent In;                  % Inputs containig trajectories of weather and setpoints
persistent ny;                  % Number of zones

if isempty(x)
   	In = evalin('base','In');                   % Evaluate input struct from workspace
    Par = evalin('base','P.Ctrl');       % Evaluate Parameters from workspace
    ny = length(Par.PI.K_i);    
    x = zeros(1,ny);
    state_coolmode = zeros(1,ny);
    state_coolmode_AHU = 0;
end

%% Compute mass flow rates and cooling/heating state for all rooms
k = find(In.Time>t,1,'first')-1;
m_sup = NaN(ny,1);
Cool_mode = NaN(ny,1);
for nz = 1:ny
    % Determine Current Setpoint
    kpre = Par.PI.T_pre(nz)*60/Par.PI.T_s;
    Tupper = min(In.T_ref_upper.Data(k,nz),In.T_ref_upper.Data(k+kpre,nz));
    Tlower = max(In.T_ref_lower.Data(k,nz),In.T_ref_lower.Data(k+kpre,nz));
    
    % Determine mode using BangBang scheme
    if T_air(nz) > Tupper
        Cool_mode(nz) = 1; 
    elseif T_air(nz) < Tlower
        Cool_mode(nz) = 0; 
    else
        Cool_mode(nz) = state_coolmode(nz);
    end
    state_coolmode(nz) = Cool_mode(nz);  
    
    % Control Error    
    if state_coolmode(nz)
        e = T_air(nz) - Tupper;
    else
      	e = Tlower - T_air(nz);
    end
    
    % Integral part
    x(nz) = x(nz) + Par.PI.K_i(nz)*e ; 
    
    % Anti Wind up
    if x(nz) > Par.m_max(nz) 
        x(nz) = Par.m_max(nz); 
    elseif x(nz)<0 
        x(nz) = 0; 
    end
    
    % Controller Output
    m_sup(nz) = Par.PI.K_p(nz)*e + x(nz);
    
    % Check Maximum
    if m_sup(nz) > Par.m_max(nz) 
        m_sup(nz) = Par.m_max(nz); 
    elseif m_sup(nz) < 0 
        m_sup(nz) = 0; 
    end
end

%% Determine Heating/Cooling Mode of AHU using bang-bang scheme
if sum(state_coolmode) >= 3 
    state_coolmode_AHU = 1;
elseif sum(state_coolmode) <= 1 
    state_coolmode_AHU = 0;
end
    
if rem(t,7*86400) < 2*86400 || rem(t,86400) < Par.T_day(1)*3600 ...
        || rem(t,86400) >= Par.T_day(2)*3600  % NightMode if Sat or Sunday or Outside of Daytime
    if state_coolmode_AHU 
        T_sup = min(Par.T_sup_night(1),In.T_out.Data(k)); 
    else
        T_sup =  Par.T_sup_night(2);
    end
else % DayTime Weekdays
    if state_coolmode_AHU 
        T_sup = Par.T_sup_day(1); 
    else
        T_sup =  Par.T_sup_day(2); 
    end
end

%% Set Massflows to 0 if States do not match
m_sup(Cool_mode ~= state_coolmode_AHU) = 0;

%% Increment counter and stop timer
t_sol = toc;
end
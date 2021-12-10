function [H_flo, H_win, varout] = ctrl_pi(t,T_air)

persistent P;                               % Parameter struct 
persistent x_h;                           	% State variable of Integral part of PI controller
persistent k;                               % Counter for current simulation step
persistent In;                              % Inputs containig trajectories of weather and setpoints
persistent kpre;                            % Number of samples of preheating
persistent H;

if isempty(x_h)
   	In = evalin('base','In');           	% Evaluate input struct from workspace
    P = evalin('base','P');                 % Evaluate Parameters from workspace
 	k = find(In.Time>t,1,'first')-1;      	% Init timestep
    x_h = 0;
    kpre = P.Ctrl.PI.T_pre*60/P.T_s;
    H = zeros(2,1);
end

%% Compute PI for heating
if mod(k-1,P.Ctrl.Std.Flo.T_s/P.T_s) == 0
    IsCoolMode = t>= P.CoolDays(1)*86400 && t<P.CoolDays(2)*86400;
    if IsCoolMode
        e = T_air - min(In.T_ref_upper.Data(k),In.T_ref_upper.Data(k+kpre));	% Control Error cooling
    else
        e = max(In.T_ref_lower.Data(k),In.T_ref_lower.Data(k+kpre)) - T_air; 	% Control Error heating
    end
    x_h = x_h + P.Ctrl.Std.Flo.K_i*e;                                                % Backward Euler
    if x_h > P.Ctrl.Std.Flo.H_max                                                    % Anti Wind up (clamping)
        x_h = P.Ctrl.Std.Flo.H_max; 
    elseif x_h<0 
        x_h = 0; 
    end
    H(1) = P.Ctrl.Std.Flo.K_p*e + x_h;                                              % Controller Output
    if H(1) > P.Ctrl.Std.Flo.H_max                                                  % Check Maximum
        H(1) = P.Ctrl.Std.Flo.H_max; 
    elseif H(1) < 0 
        H(1) = 0; 
    end
end
    
%% Assign H_win
if mod(k-1,P.Ctrl.Std.Win.T_s/P.T_s) == 0
    switch lower(P.Ctrl.Std.Win.Mode)
        % 0 is opened, 1 is closed
        case 'off'
            H(2) = 0;     
        case 'time_unocc'  
        	if In.N_pers.Data(k)==0 && IsCoolMode
                H(2) = 1; 
            else
            	H(2) = 0;
        	end
        case 'time'
            if IsCoolMode
                H(2) = 1; 
            else
            	H(2) = 0;
            end
        case 'bb'
            if T_air > P.Ctrl.Std.Win.BbHyst
                H(2) = 1;
            elseif T_air < P.Ctrl.Std.Win.BbHyst
                H(2) = 0;
            end
        otherwise 
            error('Wrong mode!');
    end
end

%% Assign Outputs
H_flo = H(1);
H_win = H(2);

%% Increment counter and stop timer
k = k + 1;
varout = NaN(3,1);
end
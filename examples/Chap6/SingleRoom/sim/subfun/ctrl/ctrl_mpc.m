function [H_flo, H_win, varout] = ctrl_mpc(t,T_air)

tic;
persistent P;                               % Parameter Structure
persistent In;                              % Input trajectories (Setpoints etc)
persistent k;                               % Current timestep
persistent opt;                             % Solver options
persistent M;                               % Model Struct
persistent SS;                              % State Spaces of Models
persistent X;                               % States of Model
persistent dQopt;                        	% optimal solution heating/cooling of previous step
persistent dHopt;                        	% optimal solution window of previous step
persistent k_start;                         % MPC start time
persistent ix_int;                          % Integrator state
Uix = 3;                                    % Input of linear subsystem that is virtual control input
Psun_ix = 2;                              	% Input of linear subsystem that corresponds to radiation

%% Init 
if isempty(P)    
    P = evalin('base','P');                 % Evaluate Parameters from workspace
    M = load(P.ModelFile);                  % load Models to structure.
    In = evalin('base','In');               % Evaluate input struct from workspace
    k = find(In.Time>t,1,'first')-1;        % Init timestep
    % Assign Prediction and control horizon 
    k_start = round((datenum(P.Ctrl.MPC.MpcActivationTime)-datenum(P.Tstart))*1440/P.T_s)+1; % First MPC Step
    dQopt = zeros(P.Ctrl.MPC.Flo.Hu_noday,1);
    dHopt = zeros(P.Ctrl.MPC.Win.Hu_noday,1);
    SolOpt = util.struct2namevaluepairs(P.Ctrl.MPC.SolOpt);
    opt = optimoptions('quadprog',SolOpt{:});
    % Get Simulation model and initial state
    type = 'modal'; minrealtol = 1e-8;
    M.RoomModFlo.updateSs('Type',type);
    SS.Room = getSs(M.RoomModFlo);
    M.Room = M.RoomModFlo;
    M.Hcv = M.HcvModFlo;

    % Est Ic such that observer converges faster
    ix_int = findIntState(SS.Room);
  	nx = length(SS.Room.A);
    Ns = 25*nx;
    try % try getting data from workspace
        D = evalin('base','SimOut_PI');
        if iscell(D)
            D = D{1};
        end
        ixs = find(D.T_air.Time>=t,1,'first');
        Y = D.T_air.Data(ixs:ixs+Ns-1);
        U = M.Room.evalInputNonlinearity([  D.In.T_out.Data(ixs:ixs+Ns-1) ...
                                            D.In.P_sun.Data(ixs:ixs+Ns-1) ...
                                            D.In.Azimuth.Data(ixs:ixs+Ns-1) ...
                                            D.In.Elevation.Data(ixs:ixs+Ns-1) ...
                                            D.H_win.Data(ixs:ixs+Ns-1) ...
                                            D.H_flo.Data(ixs:ixs+Ns-1) ...
                                            D.T_sup_flo.Data(ixs:ixs+Ns-1)],Y);
        %U(:,Uix) = [0; diff(U(:,Uix),1)];
    catch
        Y = T_air.*ones(Ns,1); 
        U = zeros(Ns,4);
   	end
    X.Room = idModels.alg.lti.estIc(Y,U,SS.Room);
    
    % Get Matrices PSI,THETA and OMEGA from Statespace object (Y = PSI*X0 + THETA*U + OMEGA*D)
    % where X0 is current state U is vectorized input trajectory and D is vectorized trajectory of known
    % future disturbances (weather predictions)
    nv = size(SS.Room.B,2);
    Hp_max = max(P.Ctrl.MPC.Hy_noday(end),P.Ctrl.MPC.Hy_day(end));
    [SS.Room.PSI, Theta] = idModels.alg.lti.calcSsPredictors(SS.Room.A,SS.Room.B,SS.Room.C,SS.Room.D,Hp_max);
 	SS.Room.THETA(:,:,1) = Theta(:,1:nv:end); % Outside temp
    SS.Room.THETA(:,:,2) = Theta(:,2:nv:end); % Solar rad
    SS.Room.THETA(:,:,3) = Theta(:,3:nv:end); % Heating/cooling
    %SS.Room.THETA(:,:,4) = Theta(:,4:nv:end); % Diff rad
    % Init Weathermodels
    if P.Ctrl.MPC.UseWeatherPredictionModel
     	% Get Statespace matrices and convert them to sparse
     	SS.Wea = getSs(M.WeaMod,[],1e-8);
        % Init State
        ix = k:M.WeaMod.Ts/(P.T_s*60):M.WeaMod.Ts/(P.T_s*60)*length(SS.Wea.A);
        X.Wea = M.WeaMod.calcIc([In.T_out.Data(ix) In.P_sun.Data(ix)]);
        fprintf('Observer for Weathermodel initilized! y_obs = \n');
        (SS.Wea.C*X.Wea)
    end 
end

%% Determine Cool-/Heat-/Shading mode
IsCoolMode = false;
if t>=P.CoolDays(1)*86400 && t<P.CoolDays(2)*86400
    IsCoolMode = true;
end
% Get Occupancy state
IsOcc = In.N_pers.Data(k)>0;
switch lower(P.Ctrl.MPC.Win.Mode)
    case 'on'
        CtrlShading = true;
    case 'on_unocc'
        CtrlShading = ~IsOcc && P.Ctrl.MPC.CtrlShading;
    otherwise
        CtrlShading = false;       
end

UseOutputConstraints = false;
if any(strcmpi(P.Ctrl.MPC.UseOutputConstraints,{'lin' 'quad'}))
    UseOutputConstraints = true;
end

%% Init Inputintegrators
if ~isfield(X,'InputIntegrator')
  	% Init state to remember u(-1)
    X.InputIntegrator = [0; 1];
    if ~CtrlShading && strcmpi(P.Ctrl.MPC.Win.Mode,'closed')
        X.InputIntegrator(2) = 0;
    end
    if CtrlShading && IsCoolMode % In Coolmode default window pos is closed.
       X.InputIntegrator(2) = 0; 
    end
end 
%% CALC CONTROL SIGNAL
% Check wether is day and get Hy an Hu
nk = P.Ctrl.MPC.Win.T_s*60/M.RoomModFlo.Ts;
t_abs = t/86400+datenum([P.Tstart(1:4) '-1-1']); % absolute time (datenum)
t_now_str = datestr(t_abs,'mm-dd HH:MM');
[IsDay,Hy,Hu_flo,Hu_win] = getHyHu(t); % get current prediction horizon and bounds 
if Hu_win*nk > Hy 
    Hu_win = ceil(Hy/nk);
end

% Predict environmental conditions and resample using linear interp
Az = In.Azimuth.Data(k+(0:Hy));
El = In.Elevation.Data(k+(0:Hy));
if P.Ctrl.MPC.UseWeatherPredictionModel 
    tq = (0:max(Hy))*P.T_s*60;  % Timepoints to be queried
    Ns = ceil(tq(end)/M.WeaMod.Ts)+1; % Number of samples to be predicted y(1) corresponds to current time
    Wea = M.WeaMod.simulate(zeros(Ns,0),X.Wea,'Ss',SS.Wea);
    Wea = resample(Wea,tq,0:M.WeaMod.Ts:(Ns-1)*M.WeaMod.Ts,'linear');
    Tout  = Wea(:,1);
    Psun = Wea(:,2);
    Psun(Psun<0) = 0;
else % Read environmental conditions from In.Data
    Tout = In.T_out.Data(k+(0:Hy));
    Psun = In.P_sun.Data(k+(0:Hy));
end
if any(In.P_sun.Data(k+(0:Hy))>.5)
    stop = 1;
end

% Predict supply Temp
if IsCoolMode
    Tsup = 18*ones(length(Tout),1);
else
    Tsup = M.Hcv.simulate([In.T_out.Data(k); Tout(2:end)]); % Supply Temp and predicted current time
end

% Optimize
if strcmpi(t_now_str,'02-06 22:00')
    stop = 1;
end
if k >= k_start   
    % Get Parameters
  	n = M.Room.InputNonlinearity(Uix).parameters{2};    % Heating Exponent
    m = M.Room.InputNonlinearity(Uix).parameters{1};    % Valve Coefficient

    % Predict Radiation on window plane by evaluating corresponding input nonlinearity
    fu = M.Room.InputNonlinearity(Psun_ix);
    Psun_win = fu.fun([Psun Az El zeros(length(El),1)],cell2mat(fu.parameters));
    
    % Predict free response EPS(0) ... EPS(Hp)
    EPS = SS.Room.PSI(1:Hy+1,:)*X.Room ... % Free movement of System
            + SS.Room.THETA(1:Hy+1,1:Hy+1,1)*Tout ... % Response to Tout
            + SS.Room.THETA(1:Hy+1,1:Hy,Uix)*ones(Hy,1)*X.InputIntegrator(1); % Response if Heating state isnt changed; 
    if ~CtrlShading 
        EPS = EPS + SS.Room.THETA(1:Hy+1,1:Hy+1,2)*Psun_win; %+ SS.Room.THETA(1:Hy+1,1:Hy+1,4)*Psun;
    else
        THETA_sun = SS.Room.THETA(1:Hy+1,1:Hy+1,2).*Psun_win'; %+ SS.Room.THETA(1:Hy+1,1:Hy+1,4).*Psun';
        EPS = EPS + THETA_sun*ones(Hy+1,1)*X.InputIntegrator(2); 
    end
    %y = M.Room.simulate([Tout Psun In.Azimuth.Data(idx) In.Elevation.Data(idx) zeros(length(idx),1) Tsup],X.Room);
    
	% get new estimate of the optimal solution (necessary for linearization of constraints)
    dQopt = [dQopt(2:end); 0];
    dHopt = [dHopt(2:end); 0];
    
	% Resample prior solution if necessary
    IsDayPrev = getHyHu(t-P.T_s);
    if IsDay(1) ~= IsDayPrev(1)
        dQopt = resample(dQopt,1:Hu_flo,1:length(dQopt)); %resample(v,xq,x,mth)
        dHopt = resample(dHopt,1:Hu_win,1:length(dHopt));
    end
       
	% Build appropriate THETA, Itri, add Inputintegration
    %QvMax = abs((Tsup - M.RoomModFlo.simulate([Tout Psun Az El (1-X.InputIntegrator(2))*ones(Hy+1,1) ones(Hy+1,1) Tsup]))^n);
    Itri_flo = tril(ones(Hu_flo));
    THETA = SS.Room.THETA(2:Hy+1,1:Hy,Uix)*tril(ones(Hy,Hu_flo));
    
    n0 = mod(nk-k+1,nk);
    Itri_win = tril(ones(Hu_win));
    if CtrlShading
        Itri_ = [zeros(n0,Hu_win)
                 kron(Itri_win,ones(nk,1))]; 
        if size(Itri_,1)>Hy
            Itri_ = Itri_(1:Hy,:);
        else
            Itri_ = [Itri_; ones(Hy-size(Itri_,1),Hu_win)];
        end
        THETA = [THETA THETA_sun(2:Hy+1,1:Hy)*Itri_];
    end

    % Get Matrices and vectors for building the optimization problem
 	[H,fT,V0] = getHf(); 

    % Get inequality constraint matrices A*U <= b 
    [A,b] = getAb();

    % Check weather heating = off and no change in shading is best thing to do
    DoOpt = true;
    u0 = [-X.InputIntegrator(1); zeros(Hu_flo-1,1)]; 
    if CtrlShading
        if IsCoolMode
            u0 = [u0; -X.InputIntegrator(2); zeros(Hu_win-1,1)]; % close
        else
            u0 = [u0; 1-X.InputIntegrator(2); zeros(Hu_win-1,1)]; % open
        end
    end
    u0 = [u0; zeros(UseOutputConstraints,1)];
    if all(A*u0-b<=0) && P.Ctrl.MPC.Wy_day==0 && P.Ctrl.MPC.Wy_noday==0
        DoOpt = false; u_opt = u0; Cost = 0; exitflag = 1; output.iterations = 0;
    end
    
    % Do opt
    if DoOpt
        if all(H(:)==0)
            [u_opt,Cost,exitflag,output] = linprog(fT,A,b,[],[],[],[],optimoptions('linprog','Display','off'));
        else
            u_opt = dQopt; 
            if CtrlShading
                u_opt = [u_opt; dHopt];
            end
            if UseOutputConstraints
                u_opt = [u_opt; 0];
            end
            [u_opt,Cost,exitflag,output] = quadprog(H,fT,A,b,[],[],[],[],u_opt,opt);
        end
        Cost = Cost + V0 + Hy*X.InputIntegrator(1); % Cost does not contain cost produced by u(-1) 
    end
    
    % Assign dQopt and dHopt
    dQopt = u_opt(1:Hu_flo);
    if CtrlShading
        dHopt = u_opt(Hu_flo+1:Hu_flo+Hu_win);
    else 
        dHopt = zeros(Hu_win,1); % leave like it is = 1
    end
    
    % Assign gamma
    if UseOutputConstraints
        gamma = u_opt(end);
    else
        gamma = NaN;
    end
    
    %% Assign output at current time k
    T_air_pred  = [T_air; THETA*u_opt(1:end-UseOutputConstraints) + EPS(2:end)];    
    Qv = Itri_flo*dQopt + X.InputIntegrator(1); 
    H_flo = (abs(Qv)./abs(Tsup(Hu_flo+1)-T_air_pred(Hu_flo+1)).^n).^(1/m);
    H_win = Itri_win*dHopt + X.InputIntegrator(2); % 1 is opend
    % Due to linearization of constraints it might happen that constraint for H_flo is violated slightly
    H_flo(H_flo>.99) = 1; H_flo(H_flo<.01) = 0;
  	H_win(H_win>.99) = 1; H_win(H_win<.01) = 0;
    % Output only first element
    H_flo = H_flo(1);
    if H_flo>0
        stop = 1;
    end
    if n0 == 0
        H_win = 1-H_win(1); % FMU treats 0 as opened and 1 as closed
    else
        H_win = 1 - X.InputIntegrator(2);
    end
else
    Cost = NaN; gamma = NaN;
    [H_flo, H_win, t_sol] = ctrl_pi(t,T_air); 
end

%% Calculate next state x_hat(k+1) = Ax_hat(k) + Bf(u(k),y(k)) + K(y(k) - y_hat(k)) (A-Priori estimate).
% Calculate actual power (since H might exceeded 1)
Qv = M.Room.InputNonlinearity(Uix).fun([H_flo Tsup(1)],cell2mat(M.Room.InputNonlinearity(Uix).parameters),T_air);

% Update state of room model
fuy = [ In.T_out.Data(k); ....
        M.Room.InputNonlinearity(Psun_ix).fun([In.P_sun.Data(k) In.Azimuth.Data(k) In.Elevation.Data(k) H_win],cell2mat(M.Room.InputNonlinearity(Psun_ix).parameters)); 
        Qv;
        %(1-H_win)*Psun(1)
        ];
Yo = SS.Room.C*X.Room + SS.Room.D*fuy;
e_obs = T_air - Yo;
X.Room = SS.Room.A*X.Room + SS.Room.B*fuy + SS.Room.K*e_obs;

% Update states of weather models
if P.Ctrl.MPC.UseWeatherPredictionModel
	Wea_o = SS.Wea.C*X.Wea;
  	X.Wea = SS.Wea.A*X.Wea + SS.Wea.K*([In.T_out.Data(k) In.P_sun.Data(k)]' - Wea_o);
end

% Update Input Integrator state
X.InputIntegrator = [Qv; 1-H_win];

% Update Counter
k = k + 1;

%% Output Info
t_sol = toc;
if ~DoOpt
    t_sol = NaN;
end
printStats()
varout = [e_obs; SS.Room.C(ix_int)*X.Room(ix_int); t_sol];


%% Func
function [H,fT,V0] = getHf()
    % Cost for Heating Power
	cT = [ones(1,Hu_flo-1) Hy-Hu_flo+1]*Itri_flo; % Weights 
	Qu = P.Ctrl.MPC.Flo.Wu_day*ones(Hu_flo,1);
	Qu(~IsDay(1:Hu_flo)) = P.Ctrl.MPC.Flo.Wu_noday; % Weight vector heating
    if IsCoolMode
       cT = -cT; 
    end
    if CtrlShading
        cT = [cT zeros(1,Hu_win)];
     	Qu_win = P.Ctrl.MPC.Win.Wu_day*ones(Hu_win,1);
        Qu_win(~IsDay(1:nk:nk*Hu_win)) = P.Ctrl.MPC.Win.Wu_noday; % Weight vector shading
        Qu = [Qu; Qu_win];
    end
        
    % Tracking error of free movement ("change nothing")
    if IsCoolMode
        E = In.T_ref_upper.Data(k+(1:Hy)) - EPS(2:Hy+1);
    else
        E = In.T_ref_lower.Data(k+(1:Hy)) - EPS(2:Hy+1);
    end
    Qy = P.Ctrl.MPC.Wy_noday*ones(Hy,1);
    Qy(IsDay(2:Hy+1)) = P.Ctrl.MPC.Wy_day;
    H = THETA'*diag(Qy)*THETA + diag(Qu);
    H=(H+H')/2;
    fT = cT - 2*E'*diag(Qy)*THETA;
    if UseOutputConstraints % Slack variable 
    	fT = [fT 0];
        H = [H zeros(size(H,1),1) 
             zeros(1,size(H,1)) 0];
        if strcmpi(P.Ctrl.MPC.UseOutputConstraints,'lin')
            fT(end) = P.Ctrl.MPC.Wconst*Hy;
        else
            H(end,end) = P.Ctrl.MPC.Wconst*Hy;
        end
    end
    if nargout > 2
        V0 = E'*diag(Qy)*E;
    end
end

function [IsDayTime,Hy,Hu_flo,Hu_win] = getHyHu(t)
    % determine day/nighttime
    if rem(t,7*86400) < 2*86400 || rem(t,86400) < P.Ctrl.MPC.T_day(1)*3600 || rem(t,86400) >= P.Ctrl.MPC.T_day(2)*3600
       	Hy = P.Ctrl.MPC.Hy_noday;
        Hu_flo = P.Ctrl.MPC.Flo.Hu_noday;
        Hu_win = P.Ctrl.MPC.Win.Hu_noday;
        IsDayTime = false;
    else
      	Hy = P.Ctrl.MPC.Hy_day;
        Hu_flo = P.Ctrl.MPC.Flo.Hu_day;
        Hu_win = P.Ctrl.MPC.Win.Hu_day;
        IsDayTime = true;
    end
    t = t + (1:Hy)'*M.Room.Ts;
    IsDayTime = [IsDayTime(1); ~(rem(t,7*86400) < 2*86400 | rem(t,86400) < P.Ctrl.MPC.T_day(1)*3600 | rem(t,86400) >= P.Ctrl.MPC.T_day(2)*3600)];
end
   
function [A,b] = getAb()
    if CtrlShading
        du0 = [dQopt; dHopt];
    else
        du0 = dQopt;
    end
    % Lower: I*U >= 0 <-> -I*U <= 0
    % Upper: I*U <= (Tsup - Y)^n = (Tsup - THETA*U - EPS)^n -> Linearize these constraints previous optimal solution du0
    L = abs(Tsup - [T_air; EPS(2:end) + THETA*du0]);
    L = [L(1:Hu_flo-1); mean(L(Hu_flo:end-1))];
    if IsCoolMode % cooling
        A = [Itri_flo zeros(Hu_flo,Hu_win*CtrlShading)];
        b = -X.InputIntegrator(1)*ones(Hu_flo,1);
        A = [A; [-Itri_flo zeros(Hu_flo,Hu_win*CtrlShading)]- n*diag(L.^(n-1))*THETA(1:Hu_flo,:)];
        b = [b; L.^n - n*diag(L.^(n-1))*THETA(1:Hu_flo,:)*du0 + X.InputIntegrator(1)*ones(Hu_flo,1)];
    else % heating
        % lower bound Q
        A = [-Itri_flo zeros(Hu_flo,Hu_win*CtrlShading)]; 
        b = X.InputIntegrator(1)*ones(Hu_flo,1);
        % upper bound Q
      	A = [A; [Itri_flo zeros(Hu_flo,Hu_win*CtrlShading)] + n*diag(L.^(n-1))*THETA(1:Hu_flo,:)];
        b = [b; L.^n + n*diag(L.^(n-1))*THETA(1:Hu_flo,:)*du0 - X.InputIntegrator(1)*ones(Hu_flo,1)];
    end
    
 	if CtrlShading
        % lower and upper bound H
       	A = [A;
             zeros(Hu_win,Hu_flo) -Itri_win;
           	 zeros(Hu_win,Hu_flo) Itri_win]; 
       	b = [b; X.InputIntegrator(2)*ones(Hu_win,1); ones(Hu_win,1)-X.InputIntegrator(2)];
	end
    
    if UseOutputConstraints % rho >= 0 -> -rho <= 0
     	% Output Bounds and slack variable
        A = [A zeros(size(A,1),1); 
             THETA -ones(Hy,1); 
             -THETA -ones(Hy,1);
             zeros(1,size(A,2)) -1];
        b = [b; In.T_ref_upper.Data((1:Hy)+k)-EPS((1:Hy)+1); EPS((1:Hy)+1)-In.T_ref_lower.Data((1:Hy)+k); 0];
    end
end

function vq = resample(v,xq,x,mth)
    if nargin <= 2 || isempty(x); x = Hu_flo(:); end
    if nargin <= 3; mth = 'previous'; end
    if strcmpi(mth,'diff')
      	vq = zeros(length(xq),size(v,2));
        vq(x,:) = v;
    else
        if isequal(x,xq)
            vq = v;
        elseif length(x) == 1
            vq = [v; v]; % hold
        else
            vq = interp1(x,v,xq(:),mth,'extrap');
        end
    end
end

function Ss = getSs(Mod,addIntIx,MinRealTol,Type,convtosparse,ObsPoles)
    if nargin > 1 && ~isempty(addIntIx) && addIntIx>0 
    	Mod.F(addIntIx) = struct('val',[1 -1],'free',[0 0],'factorized',false,'std',NaN(1,2));
    end
  	S = Mod.Ss;
    if nargin > 2 && ~isempty(MinRealTol)
        S = minreal(S,MinRealTol);
    end
    if nargin > 3 && ~isempty(Type)
        S = canon(S,Type);
    end  
    Ss.A = S.A;
    Ss.C = S.C;
    if isfield(S.InputGroup,'Measured')
        Ss.B = S.B(:,S.InputGroup.Measured);
        Ss.D = S.D(:,S.InputGroup.Measured);
    else
        Ss.B = zeros(size(Ss.A,1),0);
        Ss.D = zeros(size(Ss.C,1),0);
    end
    if nargin > 5 && ~isempty(ObsPoles)
        [Ss.K,prec,msg] = place(Ss.A',Ss.C',ObsPoles);
        Ss.K = Ss.K';
    else
        [Ss.K,Poles,Ak] = idModels.alg.lti.kalmanGain(S,'StabilityCheck',true);
    end
    
    if nargin > 4 && convtosparse
        Ss.A = sparse(Ss.A);
        Ss.B = sparse(Ss.B);
        Ss.C = sparse(Ss.C);
        Ss.D = sparse(Ss.D);
        Ss.K = sparse(Ss.K);
    end
end

function printStats()
    fprintf('Time: %s; ',t_now_str); 
    if exitflag >=0 
        msg = ['Sol. found in ' num2str(output.iterations) ' iter.'];
    else
        msg = output.message; 
    end
    fprintf('%s; Cost=%.2f; gamma=%.2f; SolTime=%.3f s.\n',msg,Cost,gamma,t_sol);
end
end

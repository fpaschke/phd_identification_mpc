function [T_s, m_s, E_obs, t_sol] = ctrl_mpc(t,T_air)

tic;
persistent P;                   % Parameter struct containing models as well
persistent X;                   % Observer states at previous time step
persistent In;                  % Inputs containig trajectories of weather and setpoint predictions
persistent Nz;                  % Number of Zones
persistent k dk;             	% Counter for current simulation step, 
persistent SS;                  % State Matrices of model and observer Gain  
persistent E;                   % Observer error
persistent u_opt;               % Saves optimal solution for hotstart
persistent opt;                 % fmincon options struct
persistent BuilMod;             % Building Model
persistent u_int;           	% State of input integrators (only used if P.Ctrl.MPC.UseDeltaU)
persistent k_start;             % Prior k_start PI Controllers will be used
global Uopt;                    % Saves optimal control input trajectories          
global Yopt;                    % Saves optimal control output trajectories
global Xobs;                    % Saves Observed state

%% Init 
if isempty(P)    
    P = evalin('base','P');                 % Evaluate Parameters from workspace
    P.Ctrl.MPC.Models = load(P.Ctrl.MPC.ModelFile);  % load Models to parameter structure.
    Nz = length(P.Ctrl.m_max);          % Set Zone number
    In = evalin('base','In');               % Evaluate input struct from workspace
    k = find(In.Time>t,1,'first')-1;        % Init timestep
    k_start = round((datenum(P.Ctrl.MPC.MpcActivationTime)-datenum(P.Tstart))*1440/P.T_s)+1; % First MPC Step
    % Get state transition matrices, Kalman gain and set initial state.
    if P.Ctrl.MPC.UseMimoModel       
        if P.Ctrl.MPC.UseSeasonalModel
            error('Not supported yet, since model hasnt been estimated!');
        else
            BuilMod = P.Ctrl.MPC.Models.RoomMod_MO;
        end
    else
        if P.Ctrl.MPC.UseSeasonalModel
            BuilMod = P.Ctrl.MPC.Models.RoomMod_SO_Seas;
        else
            BuilMod = P.Ctrl.MPC.Models.RoomMod_SO;
        end
    end
  	% Add Integrators to input if requested
 	% NOTE: If integrators are added to the inputs then observer can be
    % unstable. However this is not the problem since the only unobservable
    % states correspond to the input integrators.
    if P.Ctrl.MPC.UseDeltaU
        BuilMod = addInputIntegrators(BuilMod,6:9); % Add integrator to inputs u(6) ... u(9) (Massflows to zones)
    end
    SS.Buil = getSs(BuilMod.Ss,[],P.Ctrl.MPC.UseSparse,[]);
%     SS.BuilObs = getSs(BuilMod.Ss,[],P.Ctrl.MPC.UseSparse,[],P.T_s,'zoh');
%     T = SS.BuilObs.A\SS.Buil.A;
%     Tinv = SS.Buil.A\SS.BuilObs.A;
%     SS.BuilObs.A = Tinv*SS.BuilObs.A*T;
%     SS.BuilObs.B = Tinv*SS.BuilObs.B;
%     SS.BuilObs.K = Tinv*SS.BuilObs.K;
%     SS.BuilObs.C = SS.BuilObs.C*T;
        
    EstStates = NaN(length(SS.Buil.A),1);
    if P.Ctrl.MPC.UseDeltaU % find input integrators in statespace model.
        u_int = zeros(Nz,1); 
        ix_int = findIntState(SS.Buil,1e-8,'input');
        EstStates(ix_int) = u_int; % Input integrators are fixed to 0
    end
    Ns = 10*length(SS.Buil.A);
    dk = P.Ctrl.MPC.T_s/P.T_s;
  	try % try getting data from workspace 
        D = evalin('base','SimOut_PI');
        if iscell(D)
            D = D{1};
        end
        ixs = find(D.T_air.Time>=t,1,'first');
        ix = ixs:dk:min(ixs+dk*Ns-1,length(D.T_air.Time));
        Y = D.T_air.Data(ix,1:BuilMod.OutputDimension);
        U = BuilMod.evalInputNonlinearity([	D.In.T_out.Data(ix) ...
                                            D.In.P_sun.Data(ix) ...
                                            D.In.Azimuth.Data(ix) ...
                                            D.In.Elevation.Data(ix) ...
                                            D.T_sup.Data(ix) ...
                                            D.m_sup.Data(ix,:)],Y);
        %U(:,Uix) = [0; diff(U(:,Uix),1)];
    catch
        Y = T_air(1:BuilMod.OutputDimension)'.*ones(Ns,BuilMod.OutputDimension); 
        U = zeros(Ns,size(BuilMod.B,2));
   	end
    X.Buil = idModels.alg.lti.estIc(Y,U,SS.Buil,false,BuilMod);
    fprintf('Observer for Roommodel initilized! y_obs = \n');
    (SS.Buil.C*X.Buil)
    
    % Init Weathermodels
    if P.Ctrl.MPC.UseWeatherPredictionModel
     	% Get Statespace matrices and convert them to sparse
     	SS.Wea = getSs(P.Ctrl.MPC.Models.WeaMod.Ss,[],P.Ctrl.MPC.UseSparse,1e-8);
        % Init State
        ix = k:P.Ctrl.MPC.Models.WeaMod.Ts/(P.T_s*60):P.Ctrl.MPC.Models.WeaMod.Ts/(P.T_s*60)*length(SS.Wea.A);
        X.Wea = P.Ctrl.MPC.Models.WeaMod.calcIc([In.T_out.Data(ix) In.P_sun.Data(ix)]);
        fprintf('Observer for Weathermodel initilized! y_obs = \n');
        (SS.Wea.C*X.Wea)
    end
    % Assign Prediction and control horizon 
    if isscalar(P.Ctrl.MPC.Hu) && P.Ctrl.MPC.Hu~=0; P.Ctrl.MPC.Hu = 0:(P.Ctrl.MPC.Hu-1); end
    if isscalar(P.Ctrl.MPC.Hy); P.Ctrl.MPC.Hy = 1:P.Ctrl.MPC.Hy; end
    assert(P.Ctrl.MPC.Hu(1)==0,'Control Horizon needs to start with 0 (Current time)!'); 
    % Init optimal solution: rows: time colums: [T_sup m_1 ... m_Nz]
    [~,~,~,~,~,~,~,~,~,ix_u_day] = getBnds(t);
    u_opt = getDefaultU(ix_u_day);
    SolOpt = util.struct2namevaluepairs(P.Ctrl.MPC.SolOpt); 
    opt = optimoptions('fmincon',SolOpt{:});  
    Uopt = timeseries(NaN(Nz+1,length(P.Ctrl.MPC.Hu),(P.N_w*7*24*60)/P.Ctrl.MPC.T_s+1),(0:P.N_w*7*24*60/P.Ctrl.MPC.T_s)');
 	Uopt.TimeInfo.Units = 'hours'; 
    Uopt.TimeInfo.StartDate = P.Tstart;
    Yopt = timeseries(NaN(Nz,length(P.Ctrl.MPC.Hy),(P.N_w*7*24*60)/P.Ctrl.MPC.T_s+1),(0:P.N_w*7*24*60/P.Ctrl.MPC.T_s)');
    Yopt.TimeInfo = Uopt.TimeInfo;
    Xobs = timeseries(NaN(length(X.Buil),(P.N_w*7*24*60)/P.Ctrl.MPC.T_s+1)',(0:P.N_w*7*24*60/P.Ctrl.MPC.T_s)');
    Xobs.TimeInfo = Uopt.TimeInfo;
end

%% Optimize
isTs = mod(k-1,P.Ctrl.MPC.T_s/P.T_s) == 0;
IsDay = getBnds(t);
IsDayNext = getBnds(t+P.Ctrl.PI.T_pre(1)*3600);
if k >= k_start && ~IsDay && ~(~IsDay && IsDayNext)   % pass to PI Controller 1 Sampling step earlier
    if isTs % calculate new u_opt
        DoOpt = 1;  % Control Flag
        [~,Hy,Hu,A_bnd,b_bnd,lby,uby,lbu,ubu,ix_u_day] = getBnds(t); % get current prediction horizon and bounds

        % Predict environmental conditions and resample using linear interp
        if P.Ctrl.MPC.UseWeatherPredictionModel 
            tq = (0:max(Hy))*BuilMod.Ts;  % Timepoints to be queried
            Ns = ceil(tq(end)/P.Ctrl.MPC.Models.WeaMod.Ts)+1; % Number of samples to be predicted y(1) corresponds to current time
            Wea = P.Ctrl.MPC.Models.WeaMod.simulate(zeros(Ns,0),X.Wea,'Ss',SS.Wea);
            Wea = resample(Wea,tq,0:P.Ctrl.MPC.Models.WeaMod.Ts:(Ns-1)*P.Ctrl.MPC.Models.WeaMod.Ts,'linear');
            Tout  = Wea(:,1);
            Psun = Wea(:,2);
            Psun(Psun<0) = 0;
        else % Read environmental conditions from In.Data
            Tout = In.T_out.Data(k+dk*(0:max(Hy)));
            Psun = In.P_sun.Data(k+dk*(0:max(Hy)));
        end

        % Check weather Off is best solution 
        u0 = u_opt; 
        if P.Ctrl.MPC.UseDeltaU
            u0(1,2:end) = - u_int;
            u0(2:end,2:end) = 0;
        else
            u0(:,2:end) = 0;
        end
        C = constFun([u0(:); zeros(Nz,1)]);
        if all(C<=0); DoOpt = 0; Cost = 0; Eps = 0; u_opt = u0; end

        % Run optimization if necessary
        if DoOpt
            u_opt = [u_opt(2:end,:); u_opt(end,:)]; % Use Solution from privious step
            if P.Ctrl.MPC.UseDeltaU
                u_opt(end,2:end) = 0; 
            end
            try 
                [u_opt,Cost,exitflag,output] = fmincon(@costFun,[u_opt(:); zeros(Nz,1)],A_bnd,b_bnd,[],[],[],[],@constFun,opt);
            catch % In case something goes wrong during fmincon (eg: Nonlinear constraint function is undefined at initial point)
                exitflag = -1; 
            end
         	Eps = u_opt(end-Nz+1:end);
            if any(Eps>=0.1) || exitflag <0 % retry using default solution
                fprintf('Comfort bounds were violated or converged to infeasible point! -> Retrying Optimization using default solution \n!');
                u_opt = getDefaultU(ix_u_day);
                [u_opt,Cost,exitflag,output] = fmincon(@costFun,[u_opt(:); zeros(Nz,1)],A_bnd,b_bnd,[],[],[],[],@constFun,opt);
            end
            Eps = u_opt(end-Nz+1:end);
            u_opt = reshape(u_opt(1:end-Nz),length(Hu),Nz+1);
            ix = k+(0:max(Hy))*dk;
         	U = [Tout Psun In.Azimuth.Data(ix) In.Elevation.Data(ix) [resample(u_opt,0:Hu(end),Hu); u_opt(end,:).*ones(Hy(end)-Hu(end)-1,1); zeros(1,5)]];
         	Y = BuilMod.simulate(U,X.Buil,'Ss',SS.Buil); 
            Uopt.Data(:,:,t/(P.Ctrl.MPC.T_s*60)+1) = u_opt';
            Yopt.Data(:,:,t/(P.Ctrl.MPC.T_s*60)+1) = Y(2:end,1:Nz)';
            Xobs.Data(t/(P.Ctrl.MPC.T_s*60)+1,:) = X.Buil;
%             if k>=round((datenum('2021-1-9')-datenum(P.Tstart))*1440/P.T_s)+1
%                 lcol = lines(2);
%                 figure('Color','w','Position',[100 100 900 900]);
%                 for nz = 1:4
%                     s(nz) = subplot(5,1,nz); stairs(Y(2:end,nz),'Color',lcol(2,:)); ylabel(['$\vartheta_' num2str(nz) ' [^\circ \mathrm{C}$]']); hold on; stairs(Hy,[lby(:,nz) uby(:,nz)],'LineStyle','--','Color',lcol(2,:)); hold on; 
%                     yyaxis right; stairs(Hu,u_opt(:,nz+1),'Color',lcol(1,:)); ylabel(['$\dot{m}_' num2str(nz) '$ [kg/s]'],'Interpreter','latex'); s(nz).YColor = lcol(1,:); yyaxis left; s(nz).YColor = lcol(2,:);
%                     s(nz).XTickLabel = {}; xlim([0 96]);
%                 end
%                 s(5) = subplot(5,1,5); stairs(Hu,u_opt(:,1),'Color',lcol(2,:)); ylabel(['$\vartheta_\mathrm{Zul}$ [$^\circ \mathrm{C}$]']); xlim([0 96]); xlabel('Pr\"{a}diktionshorizont'); hold on; stairs(Hu,[lbu(:,1) ubu(:,1)],'LineStyle','--','Color',lcol(2,:)); hold on; ylim([14 31]); 
%               	yyaxis right; stairs(0:Hy(end),Tout,'Color',lcol(1,:)); ylabel(['$\vartheta_\mathrm{Aul} [^\circ \mathrm{C}$]']); s(5).YColor = lcol(1,:); yyaxis left; s(5).YColor = lcol(2,:);
%                 util.formatFigure(14); linkaxes(s,'x'); 
%                 util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\MPC\LB_MpcTrajWinter_Ts1_Hp4d.pdf','AxPosOffset',[0 -.02 0 0],'FigPosOffset',[0 0 10 -200],'SubplotYSpace',.03,'EqualX',true)
%             end
        end
    end
    %% Assign output at current time k
    T_s = u_opt(1,1);
    m_s = u_opt(1,2:end)'; 
    if P.Ctrl.MPC.UseDeltaU
        m_s = m_s + u_int;
        u_int = m_s;
    end
else
    DoOpt = 0; Cost = NaN; Eps = NaN;
    [T_s, m_s, ~, ~] = ctrl_pi(t,T_air);
    u_opt = [T_s m_s'; 22*ones(size(u_opt,1)-1,1) zeros(size(u_opt,1)-1,Nz)];
    if isTs
        stop = 1;
    end
end


%% Calculate next state x_hat(k+1) = Ax_hat(k) + Bf(u(k),y(k)) + K(y(k) - y_hat(k)) (A-Priori estimate).
t_sol = toc;
if isTs
    Ym = T_air(1:Nz+P.Ctrl.MPC.UseMimoModel); % Measurement at current timestep t
    if P.Ctrl.MPC.UseDeltaU
        Ym = [Ym; u_int]; 
    end
    fuy = BuilMod.evalInputNonlinearity([In.T_out.Data(k) In.P_sun.Data(k) In.Azimuth.Data(k) In.Elevation.Data(k) u_opt(1,:)],Ym')';
    Yo = SS.Buil.C*X.Buil + SS.Buil.D*fuy;
    X.Buil = SS.Buil.A*X.Buil + SS.Buil.B*fuy + SS.Buil.K*(Ym - Yo);
    E = Ym(1:Nz) - Yo(1:Nz);
    E_obs = E;

    if P.Ctrl.MPC.UseWeatherPredictionModel
        Wea_o = SS.Wea.C*X.Wea;
        X.Wea = SS.Wea.A*X.Wea + SS.Wea.K*([In.T_out.Data(k) In.P_sun.Data(k)]' - Wea_o);
    end
    
    %% Output Info
    printStats()
else
    E_obs = E;
end
k = k + 1;
%% Func
function [V,dV] = costFun(U)
    % Comp V
    epsilon = U(end-Nz+1:end);
    nHu = length(Hu);
    U = reshape(U(1:end-Nz),nHu,Nz+1); 
    if P.Ctrl.MPC.UseDeltaU % Go from dU to U
        U(:,2:end) = cumsum(U(:,2:end),1) + u_int';
    end
    Hp = max(Hy); 
    m_Hu = sum(U(:,2:end),2);
    m = resample(m_Hu,0:Hp-1); % resampled to 0...Hp-1
    T_sup = resample(U(:,1),0:Hp-1);
    dT = (T_sup - Tout(1:Hp)); 
    r = (1-P.Eta_heatrecovery)*P.C_air;
    Ph = r*m.*abs(dT);
    [Pel,~] = P.P_fan(m,P.Eta_fan);
    V = P.T_s*(P.Cost_el*sum(Pel) + P.Cost_heat*sum(Ph))/60 + P.Ctrl.MPC.Wconst*Hy(end)*sum(epsilon.^2);
	
    % Comp Gradient dV 
    if nargout > 1
        w = [Hu Hy(end)];
        W = diff(w)'; % Weights
        [~,dP] = P.P_fan(m_Hu,P.Eta_fan);
        dV = zeros(size(U));
        d = (dT>0) - (dT<0); % 1 if Tsup>Tout -1 else
        if P.Ctrl.MPC.UseDeltaU
            dEl = zeros(nHu,1);
            for nhu = 1:nHu
               qiT = [ones(1,nhu) zeros(1,nHu-nhu)];             
               dEl = dEl + P.T_s*P.Cost_el/60*W(nhu)*qiT'*dP(nhu);
               dV(nhu,1) = r*P.Cost_heat*P.T_s/60*m_Hu(nhu)*sum(d(w(nhu)+1:w(nhu+1)));
               dV(:,2:end) = dV(:,2:end) + r*P.Cost_heat*P.T_s/60*qiT'*sum(abs(dT(w(nhu)+1:w(nhu+1))));
            end
        else
            dEl = P.T_s*P.Cost_el/60*W.*dP;
            for nhu = 1:nHu % Ph part 
                dV(nhu,:) = r*P.Cost_heat*P.T_s/60*[m_Hu(nhu)*sum(d(w(nhu)+1:w(nhu+1))) sum(abs(dT(w(nhu)+1:w(nhu+1))))*ones(1,Nz)];
            end
        end
        dV(:,2:end) = dV(:,2:end) + dEl;
        dV = [dV(:); 2*Nz*P.Ctrl.MPC.Wconst*Hy(end)*epsilon];
    end
end

function [C,Ceq,dC,dCeq] = constFun(U) %C(u) <= 0
    epsilon = (U(end-Nz+1:end))';
	n_Hu = length(Hu);
    n_Hy = length(Hy) ;
    U = reshape(U(1:end-Nz),n_Hu,Nz+1);
    w = [Hu Hy(end)];
    Hp = max(Hy); 
    idx = k+(0:Hp)*dk;
    if P.Ctrl.MPC.UseDeltaU
        u = [Tout Psun In.Azimuth.Data(idx) In.Elevation.Data(idx) resample(U(:,1),0:Hp) resample(U(:,2:end),0:Hp,Hu(:)+1,'diff')+u_int'];
    else
        u = [Tout Psun In.Azimuth.Data(idx) In.Elevation.Data(idx) resample(U,0:Hp)];
    end
    if nargout > 2
        [y,~,dy_du] = BuilMod.simulate(u,X.Buil,'CalcInputGradients',5:9,'Ss',SS.Buil); 
    else
        y = BuilMod.simulate(u,X.Buil,'Ss',SS.Buil); % Grad wrt [dm and Tsup]
    end
    lb = lby - y(Hy+1,1:Nz) - epsilon; 
    ub = y(Hy+1,1:Nz) - uby - epsilon;
    C = [lb(:); ub(:)];
    if nargout > 2
       	dC = zeros(numel(U)+Nz,2*n_Hy*Nz); % nc x n_Hu x Nz+1
        for nz = 1:Nz
            dC(end-Nz+nz,[(1:n_Hy)+(nz-1)*n_Hy (1:n_Hy)+(Nz+nz-1)*n_Hy]) = -1; 
        end
        for nhu = 1:n_Hu
            for nz = 1:Nz
                ixC = n_Hy*(nz-1)+1:n_Hy*nz;
                dTsup = squeeze(dy_du(Hy+1,nz,1,1:Hp));
                dC(nhu,ixC) = -sum(dTsup(:,w(nhu)+1:w(nhu+1)),2);
                dC(nhu,ixC+length(C)/2) = -dC(nhu,ixC);
                for nz2 = 1:Nz
                    ixC = n_Hy*(nz2-1)+1:n_Hy*nz2;
                    dTm = squeeze(dy_du(Hy+1,nz2,nz+1,1:Hp)); % [dy_nz2(1)/dm_nz1(0) ... dy_nz2(1)/dm_nz1(Hp-1); ... dy_nz2(Hp)/dm_nz1(Hp-1)]
                    dC(nz*n_Hu+nhu,ixC) = -sum(dTm(:,w(nhu)+1:w(nhu+1)),2);
                    dC(nz*n_Hu+nhu,ixC+length(C)/2) = -dC(nz*n_Hu+nhu,ixC);
                end   
            end
        end
    end
	Ceq = [];
    dCeq = [];
    
    %% DEBUG: CHECK PREDICTIONS
    
end

function [IsDayTime,Hy,Hu,A,b,lb_y,ub_y,lb_u,ub_u,ix_tu_day] = getBnds(t)
    % determine day/nighttime
 	Hy = P.Ctrl.MPC.Hy;
  	Hu = P.Ctrl.MPC.Hu;
    if rem(t,7*86400) < 2*86400 || rem(t,86400) < P.Ctrl.T_day(1)*3600 || rem(t,86400) >= P.Ctrl.T_day(2)*3600
        IsDayTime = false;
    else
        IsDayTime = true;
    end
    if nargout > 1
        % Input bounds
        nHu = length(Hu);
        tu = t + P.Ctrl.MPC.T_s*60*Hu;
        ty = t + P.Ctrl.MPC.T_s*60*Hy;

        ix_tu_day = locGetDayIx([tu ty(ty>tu(end))]);
        ix_tu_day = [ix_tu_day(1:length(tu)-1); any(ix_tu_day(length(tu):end-1))];
               
        lb_u = zeros(nHu,Nz+1);
        lb_u(ix_tu_day,1) = P.Ctrl.T_sup_day(1);
        lb_u(~ix_tu_day,1) = min(P.Ctrl.T_sup_night(1),In.T_out.Data(k+dk*Hu(~ix_tu_day)));
        ub_u = kron([P.Ctrl.T_sup_night(2) (P.Ctrl.m_max(:))'],ones(nHu,1));
        ub_u(ix_tu_day,1) = P.Ctrl.T_sup_day(2);
        
        % Output bounds
        lb_y = In.T_ref_lower.Data(k+Hy*dk,:);
        ub_y = In.T_ref_upper.Data(k+Hy*dk,:);
        
        % Build A and b
        Nvar = numel(u_opt);
        if P.Ctrl.MPC.UseDeltaU
            % Lower Bounds
            A = -eye(nHu,Nvar);  % Supply Temp isnt integrated
            b = -lb_u(:,1);
            A = [A; zeros(Nz*nHu,nHu) -kron(eye(Nz),blkdiag(tril(ones(nHu,nHu))))];
            b = [b; reshape(-lb_u(:,2:end)+ u_int',[],1)];
            % Upper Bounds
            A = [A; eye(nHu,Nvar)];  % Supply Temp isnt integrated
            b = [b; ub_u(:,1)];
            A = [A; zeros(Nz*nHu,nHu) kron(eye(Nz),blkdiag(tril(ones(nHu,nHu))))];
            b = [b; reshape(ub_u(:,2:end) - u_int',[],1)];
        else
            % u >= u_min <-> -u <= -u_min
            A = -eye(Nvar);
            b = -lb_u(:);
            % u <= u_max
            A = [A; eye(Nvar)];
            b = [b; ub_u(:)];    
        end
        % eps >= 0
        A = [A zeros(size(A,1),Nz);
             zeros(Nz,size(A,2)) -eye(Nz)];
        b = [b; zeros(Nz,1)];
    end
    function ix_t_day = locGetDayIx(t)
        ix_t_day = true(length(t),1);
        ix_t_day(rem(t,7*86400) < 2*86400 | rem(t,86400) < P.Ctrl.T_day(1)*3600 | rem(t,86400) >= P.Ctrl.T_day(2)*3600) = false;
    end
end

function vq = resample(v,xq,x,mth)
    if nargin <= 2 || isempty(x); x = Hu(:); end
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

function Ss = getSs(S,type,convtosparse,minrealtol,Ts,mth)
	if nargin > 3 && ~isempty(minrealtol)
        if mod(minrealtol,1)~=0
            S = minreal(S,minrealtol);
        else
            S = balred(S,minrealtol);
        end
   	end
    if nargin > 4 && ~isempty(Ts)
        S = d2d(S,Ts,mth);
    end
    if nargin > 1 && ~isempty(type)
        S = canon(S,type);
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
    Ss.K = idModels.alg.lti.kalmanGain(S,'StabilityCheck',false);
    if nargin > 2 && convtosparse
        Ss.A = sparse(Ss.A);
        Ss.B = sparse(Ss.B);
        Ss.C = sparse(Ss.C);
        Ss.D = sparse(Ss.D);
        Ss.K = sparse(Ss.K);
    end
end

function printStats()
    fprintf('%s: ',datestr(t/86400+737791,'mm-dd HH:MM')); 
    if DoOpt
        if exitflag >=0 
            msg = [num2str(output.iterations) ' iter.; ' num2str(output.funcCount) ' fevals.'];
        else
            msg = output.message; 
        end
    else
        msg = 'No opt. performed!';
    end
    fprintf(['%s; Cost=%.2f; Eps=[', repmat('%.3f, ', 1, numel(Eps)-1), '%.3f]; SolTime=%.3f s.\n'],msg,Cost,Eps,t_sol);
end

function U = getDefaultU(ix_day)
    U = kron([mean([P.Sig.T_ref.T_set_lower(:,2); P.Sig.T_ref.T_set_upper(:,2)]) (P.Ctrl.m_max(:))'/2],ones(length(P.Ctrl.MPC.Hu),1));
    if nargin >= 1
        U(~ix_day,2:end) = U(~ix_day,2:end)/5; % use 20% volumeflowrates during night 
        ixDayOn = find(diff(ix_day)==1); % 1 sample prior switching to daymode
        U(ixDayOn,2:end) = ones(length(ixDayOn),1)*(P.Ctrl.m_max(:))';
   end
end
end

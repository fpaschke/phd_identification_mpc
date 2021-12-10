function e = identifyModel(obj,y,u,opt)
%Identifies NsfPolyModel using PEM Method.
%
% [obj, e] = identifyModel(obj,y,u,varargin)
% obj [NonlinearArmaxModel]: The Model object which coefficients shall be identified
% y [cell of N x ny double Matrix]: Output Data
% u [cell of N x nu double Matrix]: Input Data
% opt [opt. Structure]: Structure Specifying the options to use for estimation -> see identifyModelOptions.m
% e [cell of N x ny x 1/Hp double Matrix]: Hp-Step resp. Hp Multistep residuals

%% Init persistent Variables
persistent oldP; %to remember the last parametervector (only in case of opti is used)
persistent oldE; %to remember the last Residuals (in case of opti is used)
persistent oldDE; %to remember the last Jacobian (in case of opti is used)
persistent poly; %save polynomials as number of Iterations 
oldP = []; oldE = []; oldDE = []; poly = [];

%% Preinit
n_alpha = sum(arrayfun(@(fi) sum(fi.free),obj.InputNonlinearity)); %number of free parameters of InputNonlinearity
ny = obj.OutputDimension;   % # of outputs
nu = size(obj.B,2);         % # of Inputs to dyn. LTI Block
obj.UpdateFlag = false;     % Deactivate updating of Ss
for nout = 1:ny
    if opt.IntegrateNoise(nout) == true
        obj.PreFilter(nout).num = conv(obj.PreFilter(nout).num,[1 -1]); 
    end
end
ForcePosPoles = ~isempty(strfind(opt.ForcePolesG,'p')); 
ForceStabPoles = ~isempty(strfind(opt.ForcePolesG,'s')); 
ForcePosZeros = ~isempty(strfind(opt.ForceZerosG,'p')); 
ForceStabZeros = ~isempty(strfind(opt.ForceZerosG,'s')); 

%% Checks
assert(~ForcePosPoles || isdiag(obj.Na),'Pos. Real Poles can currently be forced to systems with diagonal A matrix only!');
if any(obj.HasOutputFeedback' & obj.Nk==0); error('Direct feedthrough (Nk=0) and outputfeedback is not supported!'); end
if any(arrayfun(@(i) any(~util.isNumber(obj.getFiParams(i))),find(obj.HasInputNonlinearity))); error('You need to supply initial values for the parameters of the InputNonlinearities!'); end
if any(opt.Hp>1) && any(obj.HasOutputFeedback); warning('The case of Nonlinear Models with outputfeedback and Prediction Horizon >1 is not theoretically validated - take care! You can still use this option and check weather the modell fits better for long prediction horizons.'); end
if any(opt.IntegrateNoise & opt.EstimateOutputOffset); warning('No Estimation of OutputOffsets for models with integration in the Noise Channel! -> Setting EstimateOutputOffset to 0!'); opt.EstimateOutputOffset = false; end
if all(opt.Hp == 1) && opt.MultiStep; opt.MultiStep = false; end %Should be faster
%assert(isdiag(obj.Na) || all(arrayfun(@(pfi) isequal(pfi.num,1) && isequal(pfi.den,1),obj.PreFilter)),'Prefilters and nondiagonal A(q) are unsupported yet!');
assert(all(obj.Nf(:) == 0),'Identification of models with F polynomial unsupported yet!');
for nr = 1:ny
    for nc = 1:ny 
        assert(~ForcePosPoles && ~ForceStabPoles || (obj.A(nr,nc).factorized || obj.Na(nr,nc) == 0),'All A Polynomials need to be factorized if ForcePolesG contains "p" (positive) or "s" (stable)!');
    end
    for nc = 1:nu
        assert(~ForcePosZeros && ~ForceStabZeros || (obj.B(nr,nc).factorized || obj.Nb(nr,nc) == -1),'All B Polynomials need to be factorized if ForceZerosG contains "p" (positive) or "s" (stable)!');
    end
  	assert(~ForcePosPoles && ~ForceStabPoles || (obj.E(nr).factorized || obj.Ne(nr) == 0),'All E Polynomials need to be factorized if ForcePolesG is s, p or sp!');
    assert((~obj.C(nr).factorized && ~obj.D(nr).factorized) || all(opt.Hp==1),'Factorized C or D and Hp>1 is not supported yet!')
end
assert(all(~[obj.E.factorized]) && all(~[obj.C.factorized]) || strcmpi(opt.StabilizationMethod,'none'),'Factorized C or E requires StabilizationMethod = none!');

%% Check length of Datasets
dmax = max(obj.getDmax);
if strcmpi(opt.Ic,'ls')
    ix_too_short = (cellfun(@length,y) - max(dmax) - max(opt.Hp) +1) <= 1;
elseif strcmpi(opt.Ic,'backcast')
    ix_too_short = (cellfun(@length,y) - max(2*dmax) - max(opt.Hp) +1) <= 1;
end
    
if any(ix_too_short)
   warning(['The Datasets ' num2str(find(ix_too_short')) ' cannot be used for identification because they are too short (If you are using Ic="backcast", then you can try Ic="ls" since it needs less samples for initialization of residual filter)! -> Removing!']);
   y = y(~ix_too_short);
   u = u(~ix_too_short);
end

%% OutputOffset
% if Outputoffset defined setup a virtual dummy model with addition virtual input
if any(opt.EstimateOutputOffset)
    oldObj = obj;
    Bv.val = NaN; Bv.free = 1; Bv.factorized = false; Bv.std = NaN;
    Bv = [obj.B repmat(Bv,ny,1)];
    I = cell2struct(cell(length(fieldnames(obj.InputNonlinearity)),1),fieldnames(obj.InputNonlinearity),1); 
    I.input_idx = obj.InputDimension+1; I.output_idx = [];
    I = [obj.InputNonlinearity; I];
    obj = idModels.NsfPolyModel(obj.Na,[obj.Nb ones(ny,1)],[obj.Nk zeros(ny,1)],obj.Nc,obj.Nd,obj.Ne,...
                                    'InputName',[obj.InputName; 'virtual_input_for_outputoffset'],...
                                    'A',obj.A,'B',Bv,'C',obj.C,'D',obj.D,'E',obj.E,...
                                    'InputNonlinearity',I,'OutputOffset',obj.OutputOffset);
	obj.UpdateFlag = oldObj.UpdateFlag;
	u = arrayfun(@(j) [u{j} ones(size(y{j},1),1)],1:length(y),'UniformOutput',0);
end

%% Set initial values 
% opt.u_norm = max(abs(cell2mat(up)));
% opt.y_norm = max(abs(cell2mat(y)));
opt.u_norm = ones(1,nu+any(opt.EstimateOutputOffset)); 
opt.y_norm = ones(1,size(y{1},2));

if ~obj.IsParameterized
    fprintf('Some of the initial parameter values of the model are NaN! -> Initializing!\n');
 	obj = obj.initPoly(y,u,opt);
end
% Stabilize if necessary
if ~strcmpi(opt.StabilizationMethod,'none')
    obj = obj.stabilize('C',opt.StabilizationFactor);
	obj = obj.stabilize('E',opt.StabilizationFactor);
end
poly = cell(length(opt.PlotZeros),1);

% Init NoiseVariance
[p0,~,~,A,b] = obj.getPvec; % Parametervector and linear constraints of Outputnonlinearities wrt free parameters
% opt_temp = opt; opt_temp.Hp = 1*ones(obj.OutputDimension,1); opt_temp.Feedback = 'raw';
% e = ts_resid(obj,y,u,opt_temp); % 1-step prediction errors
% obj.NoiseVariance = obj.calcCovariance(e,length(p0));

%% Build Constraints to stabilize C and E
hasNlCon = false;
if strcmpi(opt.StabilizationMethod,'constraint') % for Models with constrained C or E with deg = 1 we have linear constraints
    ix_y = 1;
	for nout = 1:ny
        if obj.Nc(nout) <= 2 && obj.Nc(nout)>0 
            assert(all(obj.C(nout).free(2:end)),'All C parameters need to be free if "StabilizationMethod==Constraint"!');
            if obj.Nc(nout) == 1; Acon = [-1 1]'; Ncon = 2; else; Acon = [-1 -1; 0 1; 1 -1]; Ncon = 3; end
            ixc = sum(arrayfun(@(bi) sum(bi.free),obj.A(nout,:))) + sum(arrayfun(@(bi) sum(bi.free),obj.B(nout,:))) + ix_y;
            ixc = ixc + (0:obj.Nc(nout)-1); 
            A = [A; zeros(Ncon,size(A,2))];
            A(end-Ncon+1:end,ixc) = Acon;
            b = [b; ones(Ncon,1)];
        else
            hasNlCon = true;
        end
        
      	if obj.Ne(nout) <= 2 && obj.Ne(nout)>0 
            assert(all(obj.E(nout).free(2:end)),'All E parameters need to be free if "StabilizationMethod==Constraint"!');
            if obj.Ne(nout) == 1; Acon = [-1 1]'; Ncon = 2; else; Acon = [-1 -1; 0 1; 1 -1]; Ncon = 3; end
            ixe = sum(arrayfun(@(bi) sum(bi.free),obj.A(nout,:))) + sum(arrayfun(@(bi) sum(bi.free),obj.B(nout,:))) + sum(obj.C(nout).free) + sum(obj.D(nout).free) + ix_y;
            ixe = ixe + (0:obj.Ne(nout)-1);
            A = [A; zeros(Ncon,size(A,2))];
            A(end-Ncon+1:end,ixe) = Acon;
            b = [b; ones(Ncon,1)];
     	else
            hasNlCon = true;
        end
        ix_y = ix_y + length(obj.getPvec(ny)) - n_alpha;     
	end
end

%% Build Constraints for A/E 
if ForcePosPoles || ForceStabPoles
    ix_y = 1;
    assert(isdiag(obj.Na) && all(obj.Nf(:)==0) ,'Unsupported yet!');
    for nout = 1:ny
        % A poly
        ixA = ix_y + (0:sum(obj.A(nout,nout).free)-1);
        NpA = length(ixA);
     	% E poly
     	ixE = sum(arrayfun(@(bi) sum(bi.free),obj.A(nout,:))) + sum(arrayfun(@(bi) sum(bi.free),obj.B(nout,:))) + sum(obj.C(nout).free) + sum(obj.D(nout).free) + ix_y;
       	ixE = ixE + (0:sum(obj.E(nout).free)-1); 
        NpE = length(ixE);
        % upper bounds 
        if ForceStabPoles % Poles should have mag < 1
            A_ = zeros(NpA+NpE,length(p0));
            A_(1:NpA,ixA) = eye(NpA); 
            A_(NpA+1:NpA+NpE,ixE) = eye(NpE); 
            A = [A; A_]; 
            b = [b; ones(NpA+NpE,1)];
        end
        % lower bound
        A_ = zeros(NpA+NpE,length(p0));
        A_(1:NpA,ixA) = -eye(NpA); 
        A_(NpA+1:NpA+NpE,ixE) = -eye(NpE); 
        A = [A; A_]; 
        if ForcePosPoles
            b = [b; zeros(NpA+NpE,1)];
        else
            b = [b; ones(NpA+NpE,1)];
        end
        ix_y = ix_y + length(obj.getPvec(ny)) - n_alpha; 
    end 
end

%% Build Constraints for B
if ForcePosZeros || ForceStabZeros
    assert(ny == 1,'Untested for ny>1')
    assert(isdiag(obj.Na) && all(obj.Nf(:)==0) ,'Unsupported yet!');
    ix_y = 1;
    for nout = 1:ny
        ix_y = ix_y + sum(obj.A(nout,nout).free);
        for nin = 1:nu
            ixBi =  ix_y + obj.B(nout,nin).free(1) + (0:sum(obj.B(nout,nin).free(2:end))-1);
            NpBi = length(ixBi);
            % upper bounds
            if ForceStabZeros 
                A_ = zeros(NpBi,length(p0));
                A_(1:NpBi,ixBi) = eye(NpBi); 
                A = [A; A_]; 
                b = [b; ones(NpBi,1)];
            end
            % lower bound
            A_ = zeros(NpBi,length(p0));
            A_(1:NpBi,ixBi) = -eye(NpBi); 
            A = [A; A_]; 
            if ForcePosPoles
                b = [b; zeros(NpBi,1)];
            else
                b = [b; ones(NpBi,1)];
            end
        	ix_y = ix_y + sum(obj.B(nout,nin).free);
        end
        ix_y = ix_y + length(obj.getPvec(ny)) - n_alpha; 
    end 
end

%% Transform Ap <= b to lb <= p <= ub if possible
is_ok = false(size(A,1),1); lb = []; ub = [];
if size(A,1)>0 
    for ncon = 1:size(A,1) % Check weather Ap<=b can be transformed to lp<=p<=lb
        is_ok(ncon,1) = sum(A(ncon,:)~=0)==1;
    end
    if isempty(lb) lb = -Inf*ones(length(p0),1); end 
    if isempty(ub) ub = Inf*ones(length(p0),1); end 
    for ncon = find(is_ok)'
        np = A(ncon,:)~=0;
        ai = A(ncon,np);
        if ai<0 
            lb(np) = b(ncon)/ai; 
        else
            ub(np) = b(ncon)/ai; 
        end
    end
    A(is_ok,:) = []; b(is_ok,:) = [];
end

%% Build equality constraints for ForceSymetricA
if opt.ForceSymmetricA
    % Set equal val vectors (satisfy constraints at initial point)
    for nr = 1:ny
       for nc = 1:nr-1
            obj.A(nr,nc).val = obj.A(nc,nr).val;
       end
    end
    assert(idModels.util.isDiagPoly(obj.A),'Polynomial Matrix A seems not to be diagonal!');
	
    ixProw = zeros(ny,1);
    for nr = 2:ny
        ixProw(nr) = sum(arrayfun(@(ai) sum(ai.free),obj.A(nr-1,:))) + sum(arrayfun(@(ai) sum(ai.free),obj.B(nr-1,:))) + sum(obj.C(nr-1).free) +  sum(obj.D(nr-1).free) +  sum(obj.E(nr-1).free) + sum(arrayfun(@(ai) sum(ai.free),obj.F(nr-1,:))); 
    end
    ixProw = cumsum(ixProw);
    
    ixPa = [];
    ixPb = [];
    for nr = 1:ny     
        for nc = 1:nr-1
            ixR = sum(arrayfun(@(ai) sum(ai.free),obj.A(nr,1:nc-1)));
         	ixPa = [ixPa; ixProw(nr) + ixR + (1:sum(obj.A(nr,nc).free))']; % indices of lower left diagonal
            ixR = sum(arrayfun(@(ai) sum(ai.free),obj.A(nc,1:nr-1)));
            ixPb = [ixPb; ixProw(nc) + ixR + (1:sum(obj.A(nc,nr).free))']; % indices of upper right diagonal
        end
        
    end
    Aeq = zeros(length(ixPa),size(A,2));
    for nc = 1:length(ixPa)
        Aeq(nc,ixPa(nc)) = 1;
        Aeq(nc,ixPb(nc)) = -1;
    end
else
    Aeq = zeros(0,size(A,2));
end
beq = zeros(size(Aeq,1),1);

%% Call appropriate Optimizer, compute and set Optimal pvector
obj.DoChecks = false;      % Deactivate Checks for speed
ix = regexp(opt.Solver,'_');
tb = opt.Solver(1:ix-1);
solver = opt.Solver(ix+1:end);
addNvp = util.struct2namevaluepairs(opt.Unmatched);
if strcmpi(tb,'matlab')
    if ~strcmpi(solver,'fmincon') && (size(A,1)~=0 || hasNlCon || size(Aeq,1)~=0) % Use fmincon if Ap<=b
        warning('The Solver %s does not support constraints Ax <= b; Ax == b and no nonlinear constraints -> using "fmincon"!',solver);
     	solver = 'fmincon'; 
    end 
    options = optimoptions(solver,'MaxIterations',opt.MaxIter,'FunctionTolerance',opt.FunTolAbs,'Display',opt.Display,...
            'DerivativeCheck',opt.CheckGradient,'MaxFunctionEvaluations',opt.MaxFunEval,'StepTolerance',opt.TolX,addNvp{:});
    if strcmpi(solver,'lsqnonlin')
        if ~isfield(opt.Unmatched,'Algorithm') %try to use levenberg marquardt by default.
            options.Algorithm = 'levenberg-marquardt';
        end
        options.Jacobian = opt.UseGradient;
        [p_opt,~,~] = lsqnonlin(@ts_cost_vect,p0,lb,ub,options);
    else
        options.GradObj = opt.UseGradient;
        if strcmpi(opt.StabilizationMethod,'constraint') && hasNlCon % Nonlinear constraints 
            [p_opt,~,~] = fmincon(@nlp_ts_cost,p0,A,b,Aeq,beq,lb,ub,@fmincon_ts_const,options); 
        else % Linear constraints 
            [p_opt,~,~] = fmincon(@nlp_ts_cost,p0,A,b,Aeq,beq,lb,ub,[],options); 
        end
    end
elseif strcmpi(tb,'opti')
    if strcmpi(opt.StabilizationMethod,'constraint') && hasNlCon
        warning('The Solver %s does not support nonlinear constraints -> using "ipopt"!',solver);
     	solver = 'ipopt';
    end
    % put TolX here!
    options = optiset('solver',lower(solver),'tolafun',opt.FunTolAbs,'maxiter',opt.MaxIter,'display',opt.Display,...
                        'derivCheck',opt.CheckGradient,'maxfeval',opt.MaxFunEval,addNvp{:});
	args = {};
    if any(strcmpi(solver,{'LEVMAR' 'LM_DER' 'MKLTRNLS' 'NL2SOL'})) % since NLS Solvers expect y and y_hat we need to compute ydata    
        if opt.UseGradient; args = [args {'grad' @opti_ts_grad}]; end
        if ~opt.MultiStep
            ydata = cell(sum(opt.Hp),1);
            for nout = 1:ny
                ydata{nout,1} = cellfun(@(yi) yi(dmax+opt.Hp(nout):end,nout),y,'UniformOutput',false);
            end
        else
            k = 1; ydata = cell(sum(opt.Hp),1);
            for nout = 1:ny
                for j = 1:opt.Hp(nout)
                    ydata{k,1} = cellfun(@(yi) yi(dmax+j:end,nout),y,'UniformOutput',false);
                    k = k+1;
                end
            end
        end
        ydata = cell2mat(vertcat(ydata{:})); %[y1{1}_1 y1{2}_1 ... y1{Nsets}_1 y1{1}_2 ... y2{Nsets}_2 ... y1{1}_Hp ... y2{Nsets}_Hp]
        ydata = ydata./sqrt(length(ydata));
        args = [{'fun' @opti_ts_cost 'ydata' ydata} args];
    else 
        args = [{'fun' @nlp_ts_cost} args];
        if strcmpi(opt.StabilizationMethod,'constraint')
            nc = 0; 
            if obj.Nc>2; nc = nc + obj.Nc; end
            if obj.Ne>2; nc = nc + obj.Ne; end
            if nc > 0
                args = [args {'nl' @ipopt_ts_const -Inf*ones(nc,1) ones(nc,1)}]; 
            end
        end
        if opt.UseGradient; args = [args {'grad' @nlp_ts_grad}]; end
    end
    args = [args {'x0' p0 'A' A 'b' b 'Aeq' Aeq 'beq' beq 'bounds' lb ub 'options' options}];
    Opt = opti(args{:});
    [p_opt,~,~,~] = solve(Opt);
else
    error('Choose a valid Solver!');
end
obj.setPvec(p_opt);

%% Plot Convergence of Polynomials
if ~isempty(opt.PlotZeros)
    obj.plotZeros(poly,opt.PlotZeros);
end

%% Calc Offsets
if any(opt.EstimateOutputOffset)
    A = arrayfun(@(ai) obj.getCoeff(ai),obj.A,'UniformOutput',false);
    E = arrayfun(@(ei) obj.getCoeff(ei),obj.E,'UniformOutput',false);
    ix_p_y0 = NaN(ny,1); ixp = 0; p_off = NaN(ny,1);
    for nout = 1:ny
        ix_p_y0(nout,1) = ixp + sum(arrayfun(@(ai) sum(ai.free),obj.A(nout,:))) + sum(arrayfun(@(bi) sum(bi.free),obj.B(nout,:)));
        p_off(nout,1) = opt.y_norm(nout)/opt.u_norm(end) * obj.B(nout,end).val;
        ixp = ixp + length(obj.getPvec(nout)) - n_alpha;
        AE(nout,:) = cellfun(@(ai) conv(ai,E{nout}),A(nout,logical(opt.EstimateOutputOffset)),'UniformOutput',false);
    end
    y0 = cellfun(@sum,AE)\p_off;
    p_opt = p_opt(setdiff(1:length(p_opt),ix_p_y0)); 
    obj = oldObj;
    obj.OutputOffset(logical(opt.EstimateOutputOffset)) = obj.OutputOffset(logical(opt.EstimateOutputOffset)) + y0;
    obj = obj.setPvec(p_opt);
    u = arrayfun(@(j) u{j}(:,1:end-1),1:length(y),'UniformOutput',0);
end
%% Undo normalization
for nout = 1:ny
    for nin = 1:nu
        obj.B(nout,nin).val = obj.B(nout,nin).val.*(opt.y_norm(nout)/opt.u_norm(nin));
    end
end
opt = rmfield(opt,{'y_norm' 'u_norm'});

%% Stabilize if necessary
% for no = 1:ny
%     if obj.C(nout).factorized; rC = obj.C(no).val(2:end); else; rC = roots(obj.C(no).val); end
%     if any(abs(rC)>=1)
%         warning('Predictor is unstable! (C has roots with magniture >=1)!'); 
%         if strcmpi(opt.StabilizationMethod,'hard')
%             warning('Stabilizing C!'); 
%             obj = obj.stabilize('C',opt.StabilizationFactor); 
%         end
%     end
%     if obj.E(nout).factorized; rE = obj.E(no).val(2:end); else; rE = roots(obj.E(no).val); end
%     if any(abs(rE)>=1)
%         warning('Predictor is unstable! (E has roots with magniture >=1)!'); 
%         if strcmpi(opt.StabilizationMethod,'hard')
%             warning('Stabilizing E!'); 
%             obj = obj.stabilize('E',opt.StabilizationFactor); 
%         end
%     end
% end
%% CALC RESIDUALS AND COVARIANCE MATRIX
opt_temp = opt; opt_temp.Hp = 1*ones(obj.OutputDimension,1);
e = ts_resid(obj,y,u,opt_temp); % 1-step prediction errors
obj.NoiseVariance = obj.calcCovariance(e,length(p_opt));

%% CAL COVARIANCE WRITE TO OBJ
Info = obj.Info;
Info.CovP = obj.calcCovP(y,u,opt); 
if isnumeric(Info.CovP)
    obj.setPvec(p_opt,sqrt(diag(Info.CovP)));
end
Info.Popt = p_opt;
Info.NumEstParameters = length(p_opt); 
obj.Info = Info;
obj.UpdateFlag = true; 
obj.DoChecks = true;

%% ----------------------- INLINE FUNS -------------------------------
%% General Costfun
function [E, dE] = calcResiduals(P) %[A,B_i,C,D,E,alpha_i]
    % set Pvec
    obj = obj.setPvec(P); 
    % set NoiseVariance
    if strcmpi(opt.Feedback,'sim') && any(obj.HasOutputFeedback)
        O = opt; O.Hp = 1*ones(obj.OutputDimension,1); O.Feedback = 'raw';
        e = ts_resid(obj,y,u,O); % 1-step prediction errors
        obj.NoiseVariance = obj.calcCovariance(e,length(P));
    end
    if strcmpi(opt.StabilizationMethod,'hard') %&& any(cellfun(@(r) any(r>=1),rCE))
        obj = obj.stabilize('C',opt.StabilizationFactor);
      	obj = obj.stabilize('E',opt.StabilizationFactor);
    end
    
    if nargout > 1 
        [E,dE] = ts_resid(obj,y,u,opt); % e = y - ym        
        %[E,dE] = cellfun(@(x,y) deal(x.*opt.y_norm,y.*opt.y_norm),E,dE,'UniformOutput',false); %undo normalization
    else
    	E = ts_resid(obj,y,u,opt); % e = y - ym
    end
     
    if ~isempty(opt.PlotZeros)
        for np = 1:length(opt.PlotZeros)
            poly{np} = [poly{np}; obj.(opt.PlotZeros{np}).val];
        end
    end
end

%% Vectorized Costfunction
function [e, de] = ts_cost_vect(P)
    if nargout > 1
        [E, dE] = calcResiduals(P);
    else
        E = calcResiduals(P);
    end
    
	% Remove Nans and handle multistep criterion
    if opt.MultiStep
        e = cell(sum(opt.Hp),1); 
    else
        e = cell(ny,1); 
    end
    k = 1; de = e;
    for no = 1:ny
        if opt.MultiStep
            for hp = 1:opt.Hp(no)
                e{k,1} = cellfun(@(ei) ei(dmax+hp:end,no,hp),E,'UniformOutput',false);
                if nargout > 1
                    de{k,1} = cellfun(@(dei) dei(dmax+hp:end,no,:,hp),dE,'UniformOutput',false);
                end
                k = k+1;
            end
        else     
            e{no,1} = cellfun(@(ei) ei(dmax+opt.Hp(no):end,no),E,'UniformOutput',false);
            if nargout > 1
                de{no,1} = cellfun(@(dei) dei(dmax+opt.Hp(no):end,no,:),dE,'UniformOutput',false);
            end
        end
    end
    e = cell2mat(vertcat(e{:}));  %[e1{1}_1 e1{2}_1 ... e1{Nsets}_1 e1{1}_2 ... e2{Nsets}_2 ... e1{1}_Hp ... e2{Nsets}_Hp]
    sqrtN = sqrt(size(e,1));
    e = e./sqrtN; %MMSE
    
    if nargout > 1
        de = squeeze(cell2mat(vertcat(de{:})));
        de = de./sqrtN;
    end
    
	% save results temporarily to possibly avoid recomputation of gradients (only important if opti is used)
    if strcmpi(tb,'opti')
        oldP = P;
        oldE = e; 
        if nargout > 1
            oldDE = de; 
        end
    end
end

%% Costfun, gradient and constraints for nonlinear programming implementation
function [J, dJ] = nlp_ts_cost(P)
    if opt.UseGradient
        [e, de] = ts_cost_vect(P); % saves gradient information temporarily to persistent variable to avoid recomputation
    else
        e = ts_cost_vect(P); % saves gradient information temporarily to persistent variable to avoid recomputation
    end
    J = e'*e;
    if nargout > 1
        dJ = 2*de'*e;
    end
end

% Gradient for optis nlp
function dJ = nlp_ts_grad(P)
    if all(P == oldP) % we can save some time here :-)
        e = oldE;
        de = oldDE;
    else % we need to recompute the gradient :-/
        [e, de] = ts_cost_vect(P);
    end
    %e = cell2mat(e);
 	%de = cell2mat(de);
    dJ = 2/length(e)*de'*e;
end

% Contraints to stabilize C and E (matlabs fminunc)
function [c,ceq,dc,dceq] = fmincon_ts_const(P)
    c = ipopt_ts_const(P) - 1;
    ceq = [];
   	dc = [];
    dceq = [];
end

% Contraints to stabilize C and E (opti ipopt)
function c = ipopt_ts_const(P)
    obj = obj.setPvec(P);
    c = []; 
    if obj.Nc>2
        c = [c; 
             abs(roots(obj.C.val))'];
    end
    if obj.Ne>2
        c = [c;
             abs(roots(obj.E.val))'];
    end
end

%% Costfun and Gradient for optis lsqnonlin
function yp = opti_ts_cost(P)
    if opt.UseGradient
        [e, ~] = ts_cost_vect(P); % saves gradient information temporarily to persistent variable to avoid recomputation
    else
        e = ts_cost_vect(P);
    end
    yp = ydata - e;
end

% Gradient for optis lsqnonlin
function dyp = opti_ts_grad(P)
    if all(P == oldP) % we can save some time here :-)
        de = oldDE;
    else % we need to recompute the gradient :-/
        [~, de] = ts_cost_vect(P);
    end
    dyp = -de;
end

end
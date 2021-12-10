classdef NsfPolyModel < idModels.DynamicModel
    % Represents an Polynomial Model with with nonlinear static feedback 
    %
    %     S:  A(q^-1) y[k] = q^(-nk) E^-1(q)*B(q^-1)./F(q^-1) g(u[k],y[k]) + Lden(q^-1)*C(q^-1)./Lnum(q^-1)*D(q^-1) e[k]
    %
    % where u[k] and y[k] represent the inputvector of dimension nu and ny resp. and
    % e[k] is assumed to be white noise of dimension ny with 0 mean and covariance matrix Lambda_e. 
    % Futher g(u[k],y[k]) is an general vector valued function that maps R^nu x R^ny -> R^nv, where nv is 
    % the number of virtual inputs to the LTI-part of S. If S is a purley linear system then g will be identity
    % and u = v. If g(u,y) = g(u) then S is called Hammerstein System. If g(u,y) then S has outputfeedback.
    % Notice that B./F and C./D means elementwise division, thus both matrices have the same dimensions.
    % Model will be identified with PEM Method. See identify() and identifyModel() resp..
    %
    % Further:
    %   A [ny x ny]:	Polynomial Matrix with A(0)=I.
    %   B [ny x nv]:    Polynomial Matrix. If S has outputfeedback then B needs to contain at least one delay. i. e. B(0) = 0.
    %   C [ny x ny]: 	Polynomial Matrix. C is assumed to be diagonal with C(0) = I!
    %   D [ny x ny]: 	Polynomial Matrix. D is assumed to be diagonal with D(0) = I!
    %   E [ny x ny]:	Polynomial Matrix. E is assumed to be diagonal with E(0) = I, thus E_i,i are shared by all (virt.) input TF of the corresponding output i.
    %   F [ny x nv]:    Polynomial Matrix with F(0) = I.
    %   Lden [ny x ny]: Polynomial Matrix that represents fixed part of the Noise Model -> see PreFilter property. Lden is assumed to be diagonal with Lden(0) = I!
    %   Lnum [ny x ny]: Polynomial Matrix that represents fixed part of the Noise Model -> see PreFilter property. Lnum is assumed to be diagonal with Lden(0) = I!
    % All polynomial Matrices are defined as structures of corresponding dimension with fields:
    %   val [1 x n double]: Coefficients of corresponding polynomial. [1 2 1] means P(q^-1) = 1 + 2q^-1 + 1q^-2. 
    %   free [1 x n log.]:  Coeff. indicating weather corresponding Coefficient is treated as free parameter.
    %
    % The InputNonlinearity g(u,y) is also defined as an nv x 1 structure with the following fields: 
    %   parameters [np_i x 1 double]:   parametervector p_i of g_i (i = 1...nv)
    %   free [np_i x 1 log.]:           indicating weather corresponding parameter is treated as free parameter.
   	%	input_idx [pos int vector]:     determines which inputnumbers will be fed to g_i -> g_i(u(:,input_idx),...) 
 	%  	output_idx [pos int vector]:    determines which outputnumbers will be fed to g_i -> g_i(...,y(:,input_idx))
  	% 	A [opt. Nc x np_i double]:     	Matrix of linear constraints on parameter vector p_i: The parametervector p_i of the inputnonlinearity needs to satisfy A*p_i<=b
    % 	b [opt. Nc x 1 double]:         The parametervector p of the inputnonlinearity needs to satisfy A*p_i<=b 
    % 	fun [function handle/string]:   pointing to the function implementation of  [f,df_dp,df_du,df_dy] = g_i(u,p_i,y) (i = 1...nv) where:
    %       u [N x length(input_idx) double]:   input data to g_i(u(:,input_idx),...)
    %       p_i [np_i x 1 double]:             	Parametervector of g_i   
    %       y [N x length(output_idx) double]:  output data to g_i(u(:,input_idx),p_i,y(:,output_idx))
    %   	f [N x 1 double]:                   the evalueted Function Value 
    %       df_dp [N x n_pi]:                   gradient wrt all parameters (parameters columnwise). If the corresponding
    %                                           is no optimization variable (free = false) NaN or any other value can be 
    %                                           returned from g_i in the correspondig column of df_dp. 
    %   Only necessary for computation of gradient of output wrt to input (dy_du) which can be used e.g. for MPC purposes. 
    %   See simulate.m and simulateModel.m resp.:
    %       df_du [N x nu]: gradient of g_i w.r.t all inputs i. e [dg_i/du_1 ... dg_i/du_nu]
    %       df_dy: [N x ny] gradient of g_i w.r.t all inputs i. e [dg_i/dy_1 ... dg_i/du_ny]
    %
    %   CONSTRUCTOR:
    %   obj = NsfPolyModel(Na,Nb,Nk,Nc,Nd,Ne,Nf,varargin)
    %   obj [idModels.NsfPolyModel]: Model object.
 	%   Nk [ny x nv int]: Matrix of int defining the delays in B(q^-1)
    %   Na,Nb, ... Nf [int matrix]: Matrix of int defining the orders of corresponding Polynomials.
    %   varargin [name value pairs]: Name value pairs of any valid property of NsfPolyModel, e. g. NsfPolyModel(...,'InputName','Voltage','Ts',1)
    %
    %
  	%   NOTES/REMARKS/RESTRICTIONS:
    %   - Currently unsupported options for Multi-Output-Systems: 
    %       - If A is not diagonal then Hp>1 is unsupported  
    %   - MIMO criterion is trace(Lambda), where Lambda is the empirical NoiseCovariance of e(t|t-k). det(Lambda) is unsupported currently
    %   - Multiple MISO/MIMO systems can be merged with merge(M1,M2,...,Mn) to one Model.
    %   - Currently estimation/identification of F Polynomials is unsupported
    %   - If Model has Outputfeedback and Outputoffset, the the OutputOffeset is assumed to be inside the feedbackloop!   
    %
    %   DEPENDENCIES: 
    %   - Needs Control systems toolbox for: 
    %       - conversion of polynomial description to matlab statespace object (see updateSs())
    %      	- resampling (see resample())
    %   - Needs MATLAB's Optimization-Toolbox or Free OPTI-Toolbox (see: https://www.inverseproblem.co.nz/OPTI/%) 
    %
    % EXAMPLES (for more examples see folder ...\+idModels\@NsfPolyModel\test):
    %
    % Estimate and validate simple ARMAX Model with two inputs:
    % u = [kron(rand(1,20),ones(1,50))' kron(rand(1,50),ones(1,20))']; e = .1*randn(size(u,1),1);
    % A = conv([1 -.9],[1 -.95]); C = conv([1 -.7+.4i],[1 -.7-.4i]); B1 = [0 .1]; B2 = [0 .2]; 
    % y = filter(B1,A,u(:,1)) + filter(B2,A,u(:,2)) + filter(C,A,e) + 1e2;
    % M1 = idModels.NsfPolyModel(2,1*[1 1],[1 1],2,'Ts',1,'Name','MISO ARMAX Model');
    % M1.identify(y,u,'Solver','matlab_lsqnonlin','Ic','ls','Hp',1,'MultiStep',1,'InitMethod','arx','EstimateOutputOffset',1,'IntegrateNoise',0);
    % figure('color','w'); 
    % s = subplot(2,2,1); M1.show('Boxplot',y,u,'Axes',s,'Hp',10,'ShowExpStd','--r');
    % s = subplot(2,2,2); M1.show('Scatter',y,u,'Axes',s,'Hp',10);
    % s = subplot(2,2,3); M1.show('Histogram',y,u,'Axes',s,'Hp',10);
    % s = subplot(2,2,4); M1.show('AutoCovariance',y,u,'Axes',s);
    % M1.showBode('ShowH',true)
    %
   	% Estimate an Nonlinear ARIMAX Model with g(u,y) = (u + 2*u.^2).*(3-y) and low frequent sinusoidal disturbance d inside the feedback loop:
    % rng(0); u = kron(2*rand(20,1),ones(50,1)); N = length(u); e = .00*randn(N,1); 
    % A = conv([1 -.8],[1 -.9]); B = [0 .1]; C = [1 -.85]; 
    % na = length(A)-1; nk = find(B~=0,1,'first')-1; nb = length(B)-nk; nc = length(C)-1;
    % alpha = 2;
    % F = @(u,p,y) deal((u + p(1)*u.^2).*(3-y),... %f
    %                    u.^2.*(3-y),... %df/dp
    %                    2*p(1)*u.*(3-y),... %df/du
    %                    -(u + p(1)*u.^2)); %df/dy
    % d = .1*sin(linspace(0,2*pi,N))+0.0; % Some disturbance;
    % y(1:2,1) = 0; for k = 3:N [fu,~,~,~] = F(u(k-1:-1:k-nb),alpha,y(k-1:-1:k-nb)); y(k,1) = -A(2:end)*y(k-1:-1:k-na) + B(2:end)*fu + C*e(k:-1:k-nc) + d(k); end
    % I(1).fun = F; I(1).parameters = [.5]; I(1).free = [true]; I(1).input_idx = 1;
    % M2 = idModels.NsfPolyModel(2,1,1,1,'InputNonlinearity',I);
    % M2.identify(y,u,'IntegrateNoise',1,'Solver','matlab_lsqnonlin');
    % M2.show('Prediction',y,u,'Hp',10);
    % M2.showBode('ShowH',true)
    % M2.showInputNonlinearity
    % Np = 10; ut = linspace(0,1,Np); yt = linspace(0,3,Np);
    % for k = 1:Np; for l = 1:Np; [f(k,l),~,~,~] = F(ut(k),alpha,yt(l)); [f_hat(k,l),~,~,~] = F(ut(k),M2.InputNonlinearity.parameters,yt(l)); end; end
    % figure(); surf(yt,ut,f,'FaceAlpha', 0.2, 'FaceColor', 'b','EdgeColor', 'b'); xlabel('y'); ylabel('u'); zlabel('f(u,y)'); hold on;
    % surf(yt,ut,f_hat,'FaceAlpha', 0.2, 'FaceColor', 'r','EdgeColor', 'r','LineStyle','--'); xlabel('y'); ylabel('u'); zlabel('f(u,y)'); grid on; legend({'Actual' 'Estimated'},'Location','northeast');
    
    properties
        A                   % [ny x ny struct] with fields val and free. See above.
        B                   % [ny x nv struct] with fields val and free. See above.
        C                   % [ny x 1 struct]  with fields val and free. See above.
        D                   % [ny x 1 struct]  with fields val and free. See above.
        E                   % [ny x 1 struct]  with fields val and free. See above.
        F                   % [ny x nv struct] with fields val and free. See above.
        PreFilter           % [ny x 1 struct]  with fields num and denum. Defines fixed part of Noise Model H(q^-1) = C_i(q^-1)*PreFilter(ny).denum/D_i(q^-1)*PreFilter(ny).num.
        InputNonlinearity   % [nv x 1 struct] with fields parameters and free. See above.
        Ts                  % [double scalar] indicating the sampling time of the model.
        TimeUnit            % [string] string indicting the timeunit of the model.        
    end 
    
    properties (Constant, Hidden = true)
     	ImagTol = 1e-12;    % If any of Polynomial values has imaginary part that is <=ImagTol then imaginary part will be ignored/removed.
    end
    
    properties (SetAccess = private, Hidden = true)
        DoChecks = true;    % If false then Polynomial Values and InputNonlinearity will not be checked. 
        UpdateFlag = true;  % If true Ss and F_nargout will be updated if any property changes
        F_nargout           % [nv x 1 pos int]: number of outputs that are returned of function handles of 
    end
    
    properties (SetAccess = private)
        Ss                  % [ny x nv+ny ss]: stores matlab ss object that represents the linear part of the system: will be available only if the system is fully parameterized!     
    end
    
    properties (Dependent = true)
        Nk                  % [ny x nu pos. int] fixed delays in B
        Na                  % [ny x ny pos. int] Orders of A
        Nb                  % [ny x nu pos. int] Orders of B - delay - 1
        Nc                  % [ny x 1 pos. int] Orders of C
    	Nd                  % [ny x 1 pos. int] Orders of D
        Ne              	% [ny x 1 pos. int] Orders of E
        Nf                  % [ny x nu pos. int] Orders of F
        HasOutputFeedback   % [nu x 1 logical] Indicates weather Inputchannel to the linear system has a Outputfeedback
        HasInputNonlinearity% [nu x 1 logical] Indicates weather Inputchannel to the linear system has a InputNonlinearity
        IsParameterized     % [scalar logical] Indicates weather parametervector has NaNs.
        % Methods require taking roots which can take long if polynomial
        % orders are large. These dependent properties are now simple methods
        %Type                % [ny x 1 cellstr] Model Type such as BJ,ARMA etc.
        %HasNoiseIntegration % [ny x 1 logical] Indicates weather Outputchannel has NoiseIntegration
    end
    
    methods
        %% Constructor
        function obj = NsfPolyModel(Na,varargin)
            % obj = NsfPolyModel(Na,Nb,Nk,Nc,Nd,Ne,Nf,varargin)
          	%   obj [idModels.NsfPolyModel]: Model object.
            %   Nk [ny x nv int]: Matrix of int defining the delays in B(q^-1)
            %   Na,Nb, ... Nf [int matrix]: Matrix of int defining the orders of corresponding Polynomials.
            %   varargin [name value pairs]: Name value pairs of any valid property of NsfPolyModel, e. g. M = NsfPolyModel(...,'InputName','Voltage','Ts',1)
                        
            %% get number of outputs
            ispolyorder = @(x) isnumeric(x) && ~isempty(x) && all(mod(x(:),1)==0) && all(~isinf(x(:))); 
            if ispolyorder(Na)
                ny = size(Na,1); 
            elseif ischar(Na) && all(mod(varargin{1},1)==0) % ...('armax',2,...)
                type = Na;
                n = varargin{1}(:);
                ny = length(n);
                Na = n*ones(1,ny); Nb = n; Nk = 1; Nc = zeros(1,ny); Nd = zeros(1,ny); Ne = zeros(1,ny); % Def.: arx
                switch type
                    case 'armax'
                        Nc = n;
                    case 'bj'
                        Na = zeros(ny,ny); Nc = n; Nd = n; Ne = n;
                    case 'oe'
                        Na = zeros(ny,ny); Ne = n; 
                    otherwise
                        error('Unsupported Paameterization. Consider initializing with NsfPolyModel(Na,Nb,Nc,Nd,...)');
                end
                varargin = [Nb Nk Nc Nd Ne varargin(2:end)];  %Nb,Nc,Nd,Ne,Nf
            else
                assert(~isempty(varargin),'Couldnt determine the number of outputs!');
                if ispolyorder(varargin{1})
                    ny = size(varargin{1},1);
                elseif length(varargin) > 2 && ispolyorder(varargin{3})
                    ny = size(varargin{3},1);
                elseif length(varargin) > 3 && ispolyorder(varargin{4})
                    ny = size(varargin{4},1);
                else 
                    error('Couldnt determine the number of outputs!');
                end
                if isempty(Na)
                    Na = zeros(ny);
                end
            end
            
            %% Check if setting up model from struct (e.g. NsfPolyModel(1,2,2,2,MdlStruct))
            if all(~cellfun(@ischar,varargin)) && isstruct(varargin{end}) % no name value pairs and last input is struct -> init from structure
                varargin = [varargin{1:end-1} util.struct2namevaluepairs(varargin{end})];
            end
                
            %% get number of inputs
            ix_il = find(strcmpi(varargin,'InputNonlinearity'));
            if ~isempty(varargin) && ispolyorder(varargin{1})
                Nb = varargin{1};
                % for simplified init
                if ~isempty(ix_il) 
                    if isscalar(Nb) % case Nb is scalar all B will
                        Nb = Nb*ones(ny,length(varargin{ix_il+1}));
                    elseif isvector(Nb) && length(Nb) == ny
                        Nb = Nb*ones(1,length(varargin{ix_il+1}));
                    end
                end
                nu = size(Nb,2);
                assert(size(Nb,1) == ny,'Number of rows of Nb needs to match number of outputs!');
                if length(varargin)>1 && ispolyorder(varargin{2})
                    Nk = varargin{2};
                 	if isequal(size(Nk),[1 1]) == 1 && ~isempty(ix_il)
                        Nk = Nk*ones(ny,length(varargin{ix_il+1}));
                  	end
                else
                    Nk = ones(ny,nu);
                end
              	assert(isequal(size(Nb),size(Nk)),'Nk and Nb need to have equal dimensions!');
            else
                Nb = zeros(ny,0); Nk = Nb; nu = 0;
            end
            
            %% now init Nc,Nd...
            Nc = zeros(ny,1); Nd = zeros(ny,1); Ne = zeros(ny,1); Nf = zeros(ny,nu);
            polyords = {'Nc' 'Nd' 'Ne' 'Nf'};
            ixfirstnvp = find(cellfun(@ischar,varargin),1,'first');
            if isempty(ixfirstnvp)
                IDX = length(varargin);
            else
                IDX = min(length(varargin),ixfirstnvp-1);
            end
            
            for k = 3:IDX
                if k<6 
                    exp = [ny 1]; 
                    eval([polyords{k-2} '= varargin{' num2str(k) '}(:);']);
                else
                    exp = [ny nu]; 
                    if nu > 0 
                        eval([polyords{k-2} '= varargin{' num2str(k) '};']);
                    end
                end
                assert(isequal(size(eval(polyords{k-2})),exp),[polyords{k-2} ' needs to be a vector of dimension ' mat2str(exp) '!']);
            end

            %% Now check weather InputLabels is supplied. 
            ix = find(strcmpi(varargin,'InputName'));
            if ~isempty(ix) % Determined by number of cells of Inputlabels.
                nin = length(varargin{ix+1});
            elseif ~isempty(ix_il) && isfield(varargin{ix_il+1},'input_idx') % try to determine from inputnonlinearity
                nin = max(cell2mat(arrayfun(@(il) max(il.input_idx),varargin{ix_il+1},'UniformOutput',false)));
            else
                nin = nu;   % Determine from dimension of B
            end

            %% Now make call to superclass constructor with default in out labels
         	inname = arrayfun(@(i) ['u_' num2str(i)],1:nin,'UniformOutput',0)'; 
            outname = arrayfun(@(i) ['y_' num2str(i)],1:ny,'UniformOutput',0)'; 
            obj = obj@idModels.DynamicModel('InputName',inname,'OutputName',outname);
            
            %% Set default TimeOptions
            obj.TimeUnit = 'seconds';
            obj.Ts = 1;          
                        
            %% Now initialize A, B, C, D, E, F, InputNonlinearities, IntegrateNoise, OutputOffset and PreFilter      
            str_dummy = struct('val',1,'free',false,'factorized',false,'std',NaN);
            Astr = repmat(str_dummy,ny,ny);
            Bstr = repmat(str_dummy,ny,nu); Fstr = Bstr;
            Cstr = repmat(str_dummy,ny,1); Dstr = Cstr; Estr = Cstr;
            str_dummy = struct('num',1,'den',1);
            Pf = repmat(str_dummy,ny,1);
            ixA = find(strcmp(varargin,'A'));
            if ~isempty(ixA); Astr = varargin{ixA+1}; end
         	ixB = find(strcmp(varargin,'B'));
            if ~isempty(ixB); Bstr = varargin{ixB+1}; end
            ixC = find(strcmp(varargin,'C'));
            if ~isempty(ixC); Cstr = varargin{ixC+1}; end
            ixD = find(strcmp(varargin,'D'));
            if ~isempty(ixD); Dstr = varargin{ixD+1}; end
            ixE = find(strcmp(varargin,'E'));
            if ~isempty(ixE); Estr = varargin{ixE+1}; end
            ixF = find(strcmp(varargin,'F'));
            if ~isempty(ixF); Fstr = varargin{ixF+1}; end
            ix = [ixA ixB ixC ixD ixE ixF];
            varargin([ix ix+1]) = [];
            
            % Init not supplied polynomials
            for k = 1:ny
                if isempty(ixA)
                    for l = 1:ny 
                        if k == l
                            Astr(k,l).val = [1 NaN(1,Na(k,l))];
                        else
                            Astr(k,l).val = [0 NaN(1,Na(k,l))];
                        end
                        Astr(k,l).free = [false true(1,Na(k,l))];
                        Astr(k,l).std = NaN(1,Na(k,l)+1);
                    end
                end
                for l = 1:nu
                    if isempty(ixB)
                        Bstr(k,l).val = [zeros(1,Nk(k,l)) NaN(1,Nb(k,l))];
                        Bstr(k,l).free = [false(1,Nk(k,l)) true(1,Nb(k,l))];
                        Bstr(k,l).factorized = false;
                        Bstr(k,l).std = NaN(1,length(Bstr(k,l).val));
                    end
                    if isempty(ixF)
                        Fstr(k,l).val = [1 NaN(1,Nf(k,l))];
                        Fstr(k,l).free = [false true(1,Nf(k,l))];
                        Fstr(k,l).factorized = false;
                        Fstr(k,l).std = NaN(1,length(Fstr(k,l).val));
                    end
                end
                if isempty(ixC)
                    Cstr(k,1).val = [1 NaN(1,Nc(k))];
                    Cstr(k,1).free = [false true(1,Nc(k))];
                    Cstr(k,1).std = NaN(1,Nc(k)+1);
                end
                if isempty(ixD)
                    Dstr(k,1).val = [1 NaN(1,Nd(k))];
                    Dstr(k,1).free = [false true(1,Nd(k))];
                    Dstr(k,1).std = NaN(1,Nd(k)+1);
                end
                if isempty(ixE)
                    Estr(k,1).val = [1 NaN(1,Ne(k))];
                    Estr(k,1).free = [false true(1,Ne(k))];
                    Estr(k,1).std = NaN(1,Ne(k)+1);
                end
            end
            obj.A = Astr; obj.B = Bstr; obj.C = Cstr; obj.D = Dstr; obj.E = Estr; obj.F = Fstr; obj.PreFilter = Pf;
            
            %% Set inputNonlinearity first if supplied
            ix = find(strcmpi(varargin,'InputNonlinearity'));
            if ~isempty(ix) && ~isempty(varargin{ix+1})
                obj.InputNonlinearity = varargin{ix+1};
                varargin(ix:ix+1) = [];
            else
                if ~isempty(ix); varargin(ix:ix+1) = []; end % remove empty otherwise error when setting remaining nvps
                strargs = {'fun' 'parameters' 'free' 'A' 'b' 'input_idx' 'output_idx'}; strargs{2, 1} = cell([nu 1]);
             	Fstr = struct(strargs{:});
                for k = 1:nu
                    Fstr(k).input_idx = k;
                end
                obj.InputNonlinearity = Fstr;
            end
            
            %% Now set user supplied opt. name value pairs!
            nvp = varargin(find(cellfun(@ischar,varargin),1,'first'):end);
            obj = obj.set(nvp{:});

            %% Checks
            %assert(~isempty(obj.Nb) || all(obj.Nd == 0),'Use A Polynomial to set up an AR(MA) Model!');
        end %end constructor
        
        %% Prop getters and setters    
        %A
        function A = get.A(obj)
            A = obj.A;
        end
        function set.A(obj,A)
            if (iscell(A) || (isreal(A) && isvector(A)))
                A = obj.setPolyVals(obj.A,A);
            end            
            A = obj.checkPoly(A,'A');
            obj.A = A;
            obj.updateSs();
        end
        
        %B
        function B = get.B(obj)
            B = obj.B;
        end
        function set.B(obj,B)
            if (iscell(B) || (isreal(B) && isvector(B)))
                B = obj.setPolyVals(obj.B,B);
            end   
            B = obj.checkPoly(B,'B');
            obj.B = B;
            obj.updateSs();
        end
        
        %C
        function C = get.C(obj)
            C = obj.C;
        end
        function set.C(obj,C)
            if (iscell(C) || (isreal(C) && isvector(C)))
                C = obj.setPolyVals(obj.C,C);
            end   
            C = obj.checkPoly(C,'C');
            obj.C = C;
            obj.updateSs();
        end
        
       	%D
        function D = get.D(obj)
            D = obj.D;
        end
        function set.D(obj,D)
            if (iscell(D) || (isreal(D) && isvector(D)))
                D = obj.setPolyVals(obj.D,D);
            end  
            D = obj.checkPoly(D,'D');
            obj.D = D;
            obj.updateSs();
        end
        
        %E
        function E = get.E(obj)
            E = obj.E;
        end
        function set.E(obj,E)
            if (iscell(E) || (isreal(E) && isvector(E))) %% init from cell -> assumes factorized = 0
                E = obj.setPolyVals(obj.E,E);
            end  
            E = obj.checkPoly(E,'E');
            nu = size(obj.Nb,2);
            for ny = 1:obj.OutputDimension
                if nu == 0 && length(E(ny).val)>1
                    warning('Youve specified an model with no inputs to the LTI Block but Ne>0, which makes no sense! -> Setting E=1!')
                    E(ny).val = 1; E(ny).std = 1; E(ny).free = 0; E(ny).factorized = 0;
                end
            end
            obj.E = E;
            obj.updateSs();
        end
        
       	%F
        function F = get.F(obj)
            F = obj.F;
        end
        function set.F(obj,F)
            if (iscell(F) || (isreal(F) && isvector(F))) 
                F = obj.setPolyVals(obj.F,F);
            end  
            F = obj.checkPoly(F,'F');
            obj.F = F;
            obj.updateSs();
        end
        
        %InputNonlinearity
        function F = get.InputNonlinearity(obj)
            F = obj.InputNonlinearity;
        end
        function set.InputNonlinearity(obj,F)
            F = F(:); 
            if obj.DoChecks
                expFields = {'fun' 'parameters' 'free' 'A' 'b' 'input_idx' 'output_idx' 'std'};
                assert(length(F)==size(obj.Nb,2),'InputNonlinearity shall have the same size as inputs to the linear block (size(Nb,2))!');
                for ni = 1:length(F)
                    % set outputnonlin if not defined and convert 2 handle
                    if isempty(F(ni).fun); F(ni).fun = []; end
                    if ischar(F(ni).fun) && ~isempty(F(ni).fun); F(ni).fun = str2func(F(ni).fun); end % convert to function handle if is a string

                    % add fields if not available yet
                    if ~isfield(F(ni),'parameters'); F(ni).parameters = {}; end
                    if isnumeric(F(ni).parameters); F(ni).parameters = num2cell(F(ni).parameters); end
                    npar = sum(arrayfun(@(j) isnumeric(F(ni).parameters{j}),1:length(F(ni).parameters)));

                    if ~isfield(F(ni),'A'); F(ni).A = []; end
                    if isempty(F(ni).A); F(ni).A = zeros(0,npar); end
                    if ~isfield(F(ni),'b'); F(ni).b = []; end

                    if ~isfield(F(ni),'free') || isempty(F(ni).free); F(ni).free = true(1,npar); end      
                    if ~isfield(F(ni),'std') || isempty(F(ni).std); F(ni).std = NaN(1,npar); end

                    if ~isfield(F(ni),'output_idx'); F(ni).output_idx = []; end                    
                    if ~isempty(F(ni).fun) && nargin(F(ni).fun)>=3 && isempty(F(ni).output_idx)
                        % warning('The InputNonlinearity %i (%s) has at least 3 inputs but you havent specified output_idx! Assuming no Outputfeedback!',ni,func2str(F(ni).fun));
                        % F(ni).output_idx = 1; -> Shoundnt be done: If
                        % f(u,p,y) is used but y is only used if output_idx
                        % is defined (see eg idModels.func.fun_boiler_eff)
                    end

                    % Make coversions
                    assert(isempty(F(ni).free) || all(F(ni).free==0 | F(ni).free==1),'"free" needs to be a logic vector');
                    F(ni).free = logical(F(ni).free(:))';
                    F(ni).parameters = (F(ni).parameters(:))'; 
                    if isempty(F(ni).fun) && isempty(F(ni).input_idx) %If fun is not specified then we will assume direct feedthrough of input ni
                        F(ni).input_idx = ni;
                    end
                    if ~isempty(F(ni).fun) && (isempty(F(ni).output_idx) && nargin(F(ni).fun)>=3)
                    	warning('The InputNonlinearity %i (%s) has at least 3 inputs but you havent specified output_idx! If youre trying to setup an InputNonlinearity with outputfeedback, then make sure that you specify "output_idx"! Your current setup corresponds to no outputfeedback (Hammerstein-Inputnonlinearity)!',ni,func2str(F(ni).fun));
                    end

                    % Checks
                    assert(all(cellfun(@(fn) any(strcmp(expFields,fn)),fieldnames(F(ni)))),...
                        'All structures shall have the fields: "fun" "parameters" "free" "std" "A" "b" and "input_idx" and "output_idx"');
                    assert(isempty(F(ni).fun) || isa(F(ni).fun,'function_handle'),'"fun" shall be a function handle or empty');
                    assert(all(cellfun(@(x) isnumeric(x) && isscalar(x) || ischar(x),F(ni).parameters)),'"parameters" shall be a cell with numerical values or a string that contains a reference to a parameter!');
                    assert(npar==length(F(ni).free) && npar==length(F(ni).std),'"free" shall have same length as there are numerical values in "parameters"');
                    assert(size(F(ni).A,1) == size(F(ni).b,1),'"A" and "b" need the same number of rows!');
                    assert(((size(F(ni).A,2) == length(F(ni).free)) || isempty(F(ni).A)) && (size(F(ni).b,2)==1 || isempty(F(ni).b)),'"A" and "b" shall have the correct dimensions!');
                    assert(all(mod(F(ni).input_idx,1)==0) && all(F(ni).input_idx<=obj.InputDimension) && all(F(ni).input_idx>=1),'All "input_idx" need to pos. integers wich need to be <= then the number of inputs!');
                    assert(~isempty(F(ni).fun) || length(F(ni).input_idx) == 1,'If no inputnonlinearity is defined then input_idx must be of length 1');
                    for np = 1:npar
                        if ischar(F(ni).parameters{np})
                            ix_ = strfind(F(ni).parameters{np},'_');
                            ixf = str2double(F(ni).parameters{np}(1:ix_));
                            assert(ixf~=ni,'References to own parametervorctor do not make sense!');
                        end
                    end
                end
            end
            obj.InputNonlinearity = F;           
            obj.updateF_nargout();
        end
        
        %Ts
        function Ts = get.Ts(obj)
            Ts = obj.Ts;
        end
        function set.Ts(obj,ts)
            assert(isreal(ts) && isscalar(ts) && ts>0 && ~isinf(ts),'Ts needs to be a pos real scalar double value!');
            obj.Ts = ts;
            obj.updateSs();
        end
        
     	%TimeUnit
        function Ts = get.TimeUnit(obj)
            Ts = obj.TimeUnit;
        end
        function set.TimeUnit(obj,tu)
            tu = lower(tu);
            assert(any(strcmp(idModels.util.getSupportedTimeUnits,tu)),'The Timeunit needs to be a supported timeunit!');
            obj.TimeUnit = tu;
            obj.updateSs();
        end
        
        %PreFilter
        function Pf = get.PreFilter(obj)
            Pf = obj.PreFilter;
        end
        function set.PreFilter(obj,Pf)  
            if iscell(Pf)
%                 if obj.OutputDimension == 1 && length(Pf)==2 % 
%                     Pf = {Pf};
%                 end
                Pf = idModels.util.polyCell2struct(Pf);
            end
            if isscalar(Pf) && obj.OutputDimension>1
                Pf = repmat(Pf,obj.OutputDimension,1);
            end
            assert(all(arrayfun(@(pfi) isstruct(pfi) && all(isfield(pfi,{'num' 'den'})) && pfi.num(1) == 1 && pfi.den(1) == 1 && ...
                        all(isreal(pfi.num) & ~isnan(pfi.num)) && all(isreal(pfi.den) & ~isnan(pfi.den)),Pf)) && length(Pf) == obj.OutputDimension,...
                        'PreFilter to be set needs to be a ny x 1 struct with monic double vectors num and den!')
            obj.PreFilter = Pf(:);
            obj.updateSs();
        end
        
        %UpdateFlag
        function set.UpdateFlag(obj,flag)
            assert(isscalar(flag) && (flag==1 || flag==0),'UpdateFlag needs to be a logical!');
            ufState = obj.UpdateFlag;
            obj.UpdateFlag = logical(flag);
            if ufState == false && obj.UpdateFlag == true
            	obj.updateSs(); % Force Update
            end
        end
     	function flag = get.UpdateFlag(obj)
            flag = obj.UpdateFlag;
        end
            
        %DoChecks
        function set.DoChecks(obj,flag)
            assert(isscalar(flag) && (flag==1 || flag==0),'UpdateFlag needs to be a logical!');
            obj.DoChecks = logical(flag);
        end
     	function flag = get.DoChecks(obj)
            flag = obj.DoChecks;
        end

        %% Prop getters for dependent variables
        % na
        function na = get.Na(obj)
            na = obj.getPolyDegree(obj.A);
        end
        
        % nb
        function nb = get.Nb(obj)
            nb = obj.getPolyDegree(obj.B) - obj.Nk;
        end
        
        % nk
        function nk = get.Nk(obj)
            nk = obj.getPolyDelay(obj.B);
        end
        
        % nc
        function nc = get.Nc(obj)
            nc = obj.getPolyDegree(obj.C);
        end
        
    	% nd
        function nd = get.Nd(obj)
            nd = obj.getPolyDegree(obj.D);
        end
        
        % ne
        function ne = get.Ne(obj)
            ne = obj.getPolyDegree(obj.E);
        end
        
        % nf
        function nf = get.Nf(obj)
            nf = obj.getPolyDegree(obj.F);
        end
        
      	% HasOutputFeedback
        function hof = get.HasOutputFeedback(obj)
            hof = false(length(obj.InputNonlinearity),1);
            for ni = 1:length(obj.InputNonlinearity)
                if ~isempty(obj.InputNonlinearity(ni).output_idx); hof(ni,1) = true; end
            end
        end
        
    	% HasInputNonlinearity
        function hinl = get.HasInputNonlinearity(obj)
            hinl = logical(arrayfun(@(i) ~isempty(i.fun),obj.InputNonlinearity));
        end
        
        % IsParameterized
        function isi = get.IsParameterized(obj)
            if ~isstruct(obj.A) || ~isstruct(obj.B) || ~isstruct(obj.C) || ~isstruct(obj.D) || ~isstruct(obj.E) || ~isstruct(obj.F) % Not initialized yet: B and F can be empty if autonomous system
                isi = false;
            else
                [~,~,P] = obj.getPvec; 
                if any(isnan(P))
                    isi = false;
                else
                    isi = true;
                end
            end
        end
        
        % HasNoiseIntegration
       	function h = HasNoiseIntegration(obj,outs)
            if nargin == 1
                outs = 1:obj.OutputDimension;
            end
            tol = 1e-8;
            h = false(length(outs),1);
            for no = outs
                rD = roots(obj.PreFilter(no).num);
            	if any(abs(rD-1) <= tol) 
                    h(no,1) = true; return;
            	end 
                h(no,1) = obj.polyHasZero(obj.D(no),1);
            end
        end
        
     	% Type
        function type = Type(obj)
            type = cell(obj.OutputDimension,1);
            for no = 1:obj.OutputDimension
                if (any(obj.Nf(no,:)~=0) || any(obj.Ne(no)~=0)) && obj.InputDimension > 0
                    if obj.Nc(no) == 0 && (obj.Nd(no) == 0 || obj.Nd(no) == 1 && obj.polyHasZero(obj.D(no),1)) 
                        type{no} = 'OE';
                    else
                        type{no} = 'BJ'; 
                    end                   
                    if obj.HasNoiseIntegration(no); type{no} = [type{no} ' with integrated noise']; end
                else
                    type{no} = [];
                    if any(obj.Na(no,:)~=0); type{no} = [type{no} 'AR']; end %A at least one Ai is not 0
                    if obj.Nd(no)-obj.polyHasZero(obj.D(no),1) > 0; type{no} = [type{no} 'AR']; end %D at least one D is not 0 and not Integrated 
                    if obj.HasNoiseIntegration(no); type{no} = [type{no} 'I']; end
                    if obj.Nc(no)>0; type{no} = [type{no} 'MA']; end
                    if obj.InputDimension > 0; type{no} = [type{no} 'X']; end
                end
                if any(obj.HasOutputFeedback)
                    type{no} = ['Nonlinear ' type{no}];
                elseif any(obj.HasInputNonlinearity)
                    type{no} = ['Hammerstein ' type{no}];
                end
            end
        end
        
        %Ss
        % Ss will be updated if any of the polynomials A,B,...,F or NoiseVariance, Ts, PreFilter or TimeUnit Changes
        function ss = get.Ss(obj) 
            ss = obj.Ss;
        end
        
        % Will be updated if InputNonlinearity Changes
        function no = get.F_nargout(obj)
            no = obj.F_nargout;
        end
               
        %% External Method implementations
        printParameters(obj,varargin)
        defactorize(obj,poly,ix)
    	factorize(obj,poly,ix)
        resample(obj,Ts,mth)
        initParameters(M1,M2,varargin)
        [G_num,G_den,H_num,H_den] = calcGH(obj,varargin)
        x0 = calcIc(obj,y,u,varargin)
        [y,x,dy_du] = simulateModel(obj,u,varargin)
      	yp = calcPredictions(obj,y,u,varargin)
        M = merge(m1,varargin)
        [x_hat,y_hat,e] = observe(obj,y,u,x0)
        e = identifyModel(obj,y,u,opt)
        opt = identifyModelOptions(obj,varargin)
        f = showBode(obj,varargin)
        ax = showInputNonlinearity(obj,ix_f,u0,varargin)
    end
    
    %% public hidden
    methods (Hidden = true)
        %% Update Statespace Object
        function updateSs(obj,varargin) % Needs to be public because it is accessed from idModels.set.NoiseVariance
        % Builds a StateSpace object for the linear block of the system.
        % 
        % Ss = getSs(obj,type)
        % type [opt. string]: 'modal' or 'companion'. See canon()
        % Ss [ss]: Matlab statespace object.
            if obj.UpdateFlag && obj.IsParameterized && all((~isnan(obj.NoiseVariance(:)))) && ~isempty(obj.NoiseVariance)                
                p = inputParser(); 
                addParameter(p,'MinRealTol',[],@(x) isempty(x) || isscalar(x) && x>=0)
                addParameter(p,'Type',[],@(x) isempty(x) || any(strcmpi(x,{'modal' 'companion'})));
                addParameter(p,'BalRedOrder',[]);
                addParameter(p,'NumZeroTol',[],@(x) isempty(x) || isscalar(x) && x>=0);
                addParameter(p,'Sparse',false,@(x) x == 1 || x == 0);
                parse(p,varargin{:});
                
                nu = size(obj.Nb,2);
                [G_num, G_den, H_num, H_den] = obj.calcGH('AbsorbNoiseVariance',true);
                GH = tf([G_num H_num],[G_den H_den],obj.Ts,'Variable','z^-1','TimeUnit',obj.TimeUnit);
                Ssobj = ss(GH);
                Ssobj.InputGroup.Measured = 1:nu;
                Ssobj.InputGroup.Noise = nu+1:nu+obj.OutputDimension;
                % Ss = ss(idpoly(arrayfun(@(x) x.val,obj.A,'UniformOutput',0),arrayfun(@(x) x.val,obj.B,'UniformOutput',0),...
                %         arrayfun(@(x) x.val,obj.C,'UniformOutput',0),arrayfun(@(x) x.val,obj.D,'UniformOutput',0),[],obj.NoiseVariance,obj.Ts,'TimeUnit',obj.TimeUnit),'augmented');
                Ssobj = absorbDelay(Ssobj);
                if ~isempty(p.Results.MinRealTol)
                    Ssobj = minreal(Ssobj,p.Results.MinRealTol); 
                end
                if ~isempty(p.Results.BalRedOrder)
                    Ssobj = balred(Ssobj,p.Results.BalRedOrder); 
                end
%                 if ~isempty(p.Results.NumZeroTol)
%                   	Ssobj.A(abs(Ssobj.A)<p.Results.NumZeroTol) = 0;
%                     Ssobj.B(abs(Ssobj.B)<p.Results.NumZeroTol) = 0;
%                     Ssobj.C(abs(Ssobj.C)<p.Results.NumZeroTol) = 0;
%                     Ssobj.D(abs(Ssobj.D)<p.Results.NumZeroTol) = 0;
%                 end
            	if ~isempty(p.Results.Type)
                    Ssobj = canon(Ssobj,p.Results.Type);
                end
                if p.Results.Sparse
                	Ssobj = sparss(Ssobj);
                end
             	obj.Ss = Ssobj;
            end
        end
        
        %% Determines number of outputs of g(u,y)
        function updateF_nargout(obj)
            % Find out how many outputs will be returned by f. Nargout might fail in cases where the inputnonlinearity
            % is defined by the deal function! Thus this messy workaround is needed...
          	if obj.UpdateFlag
                f = obj.InputNonlinearity;
                nu = length(f);
                nouts = NaN(nu,1);
                for ni = 1:nu
                    fi = f(ni);
                    if isempty(fi.fun)
                        nouts(ni) = 1;
                    else
                        nouts(ni) = nargout(fi.fun);
                    end
                    if nouts(ni) == -1 % could be deal function
                        no = 0;
                        while 1
                            outs = cell(1,no);
                            try
                                if obj.HasOutputFeedback(ni)
                                    [outs{:}] = fi.fun(zeros(1,length(fi.input_idx)),obj.getFiParams(ni),zeros(1,length(fi.output_idx)));
                                elseif obj.HasInputNonlinearity(ni)
                                    [outs{:}] = fi.fun(zeros(1,length(fi.input_idx)),obj.getFiParams(ni));
                                else
                                    no = 1;
                                end
                                nouts(ni) = no;
                                break;
                            catch
                                no = no + 1;
                                if no > 10
                                   error('InputNonlinearity cannot be evaluated! Check the definition!'); 
                                end
                            end
                        end
                    end
                end
                obj.F_nargout = nouts;
            end
        end    
        
        %% Calc k-step Predictor
        function [F,qkG,dF_da,dF_dc,dF_dd] = calcKStepPredictor(obj,k) %called in show thus public
            % F [ny x ny cell]: 
            % dF_da: [dF11_dA11 0 ...
            %         0 dF21_dA22 ...] % Assumes diagonal A
            % dF_dc: [dF11_dC1 dF12_dC2 ...
            %         dF21_dC1 dF21_dC2 ...] % Assumes diagonal C
            % dF_dd: [dF11_dD1 dF12_dD2 ...
            %         dF21_dD1 dF21_dD2 ...] % Assumes diagonal D

            ny = obj.OutputDimension;
            if isscalar(k)
               k = ones(1,obj.OutputDimension)*k; 
            end
            AA = idModels.util.polyStruct2cell(obj.A,'DefactorizePoly',true);
            if isdiag(obj.Na) || all(k==1) 
                Ai_num = num2cell(eye(obj.OutputDimension));
                Ai_den = AA; 
                for no1 = 1:ny; for no2 = 1:ny; if no1~=no2; Ai_den{no1,no2} = 1; end; end; end
            else
                error('Unsupported yet!');
                % Ga = tf(idpoly(AA,[],[],[],[],1),'Noise'); % Compute A^-1
                % Ai_num = Ga.num;
                % Ai_den = Ga.den;
            end
            na_den = cellfun(@length,Ai_den)-1;

            F = cell(ny,ny); qkG = F; dF_da = F; dF_dc = F; dF_dd = F;
            % H = A^-1 C*Np./ D*Zp
            for r = 1:ny % row of H
                for c = 1:ny % column of H
                    %% CALC F and G
                    n_pf_den = length(obj.PreFilter(c).den)-1;
                    n_pf_num = length(obj.PreFilter(c).num)-1;
                    
                    n = max(na_den(r,c)+obj.Nd(c)+n_pf_num,obj.Nc(c)+na_den(c)+n_pf_den); % maximum zähler und nennergrad
                    ng = max(na_den(r,c)+obj.Nd(c)-1+n_pf_num,obj.Nc(c)+n_pf_den-k(r));
                    ADZ = conv(conv(Ai_den{r,c},obj.D(c).val),obj.PreFilter(c).num); 
                    A_ = [ADZ zeros(1,n-na_den(r,c)-obj.Nd(c)-n_pf_num)]; 
                    CN = conv(conv(Ai_num{r,c},obj.C(c).val),obj.PreFilter(c).den); ncp = length(CN) - 1; 
                    C_ = [CN zeros(1,n-ncp)]; 
                    [F{r,c},qkG_] = deconv([C_ zeros(1,k(r)-1)],A_);
                    qkG{r,c} = qkG_(1:end-(n-ng-1)); %remove zeros

                    %% CALC Deriv
                    if nargout > 2 % dF/dai
                        W = conv(obj.D(c).val,obj.PreFilter(c).num); nw = length(W)-1;
                        N = conv(W,Ai_den{r,c}); nn = length(N)-1;
                        dF_da{r,c} = zeros(obj.Na(r,c),k(r)); % dF/dai(nr,nc): parameter zeilenweise
                        if obj.A(r,c).factorized; [~,dA] = idModels.util.defactorizePoly(obj.A(r,c).val); end
                        for l = 1:obj.Na(r,c) %diff Index
                            for j = 1:k(r)-1 %length(F)-1 F index
                                x = 0; 
                                for i = 1:min(j,n) 
                                    if obj.A(r,c).factorized
                                        wi = 0;
                                        for k_ = 0:min(i,obj.Na(r,c))
                                            if i-k_>=0 && i-k_<=nw; wii = W(i-k_+1); else; wii = 0; end
                                            wi = wi + wii*dA(l+1,k_+1);
                                        end
                                    else
                                        if i-l>=0 && i-l<=nw; wi = W(i-l+1);  else; wi = 0; end
                                    end
                                    if i<=nn; ni = N(i+1); else; ni = 0; end
                                    x = x + ni*dF_da{r,c}(l,j-i+1) + wi*F{r,c}(j-i+1);
                                end
                                dF_da{r,c}(l,j+1) = -x;
                            end
                        end
                    end

                    if nargout > 3 % dF/dci(nc)
                        Np = conv(obj.PreFilter(c).den,Ai_num{r,c}); nnp = length(Np)-1;
                        dF_dc{r,c} = zeros(obj.Nc(c),k(r)); % parameter zeilenweise
                        for l = 1:obj.Nc(c)
                            for j = 1:k(r)-1
                                x = 0;
                                for i = 1:min(j,n)
                                    if i<=nn; ni = N(i+1); else; ni = 0; end
                                    % if i<=nnp; npi = Np(i+1); else; npi = 0; end
                                    x = x + ni*dF_dc{r,c}(l,j-i+1); 
                                end
                                if j>=l && j-l<=nnp
                                    dF_dc{r,c}(l,j+1) = Np(j-l+1) - x;
                                else
                                    dF_dc{r,c}(l,j+1) = - x;
                                end
                            end
                        end
                    end

                    if nargout > 4 % dF/ddi(nc)
                        W = conv(Ai_den{r,c},obj.PreFilter(c).num); nw = length(W)-1;
                        dF_dd{r,c} = zeros(obj.Nd(c),k(r)); % parameter zeilenweise
                        for l = 1:obj.Nd(c)
                            for j = 1:k(r)-1 %length(F)-1
                                x = 0; 
                                for i = 1:min(j,n) 
                                    if i-l>=0 && i-l<=nw; wi = W(i-l+1);  else; wi = 0; end
                                    if i<=nn; ni = N(i+1); else; ni = 0; end
                                    x = x + ni*dF_dd{r,c}(l,j-i+1) + wi*F{r,c}(j-i+1);
                                end
                                dF_dd{r,c}(l,j+1) = -x;
                            end
                        end
                    end
                end
            end
        end
        
        %% Calculate Parameter Covariance
        function covPk = calcCovP(obj,y,u,opt)
            %CALCCOVP Calculate Covariance Matrix of parameter estimates.
            %
            % covPk = calcCovP(obj,y,u)
            % covPk [np x np double]: Empirical covariance of parameter estimates.
            % obj [idModels.NsfPolyModel]: Instance of object.
            % y [cell of N x ny double]: OutputData
            % u [cell of N x nu double]: InputData

            %% Init
            Hp = opt.Hp(:);
            dmax = max(obj.getDmax-Hp+1);
            lambda = obj.NoiseVariance;
            ny = obj.OutputDimension;
            if ~iscell(y); y = {y}; end
            if ~iscell(u); u = {u}; end
            if (ny > 1 && any(Hp>1))  
                covPk = 'Couldnt be estimated because k- and multi-step criterion isnt supported for MIMO systems yet!';
                return;
            elseif opt.MultiStep == true && any(Hp > 1)
                % warning('Function does not support Multistep Criterion and/or MIMO case yet!')
                covPk = 'Couldnt be estimated because multistep criterion isnt supported yet!';
                return;
            end
            Nsets = length(y);
            Nss = cellfun(@length,y);            
            
            %% Calculate Gradients
            [~,de] = obj.ts_resid(y,u,opt);
            npar = size(de{1},3);

            %% Calc E(de(t),de(t-j)): Matrices will be avarged wrt to the length of the dataset
            if all(Hp == 1) && ny > 1 % assumes trace criterion cost function here!
                for ns = 1:Nsets % Remove NaNs
                    de{ns} = de{ns}(dmax+1:end,:,:);
                end
              	Psi = cell2mat(de); 
                if verLessThan('matlab','9.9') % for matlab version less then 2020b pagewise multiplication needs to be done manually
                    N = size(Psi,1);
                    Q = NaN(N,npar,npar);
                    R = Q;
                    for ns = 1:N % %Pagewise array multiplication
                        H = shiftdim(Psi(ns,:,:),1);
                        Q(ns,:,:) = H'*H;
                        R(ns,:,:) = H'*lambda*H;
                    end     
                    Q = squeeze(mean(Q,1));
                    R = squeeze(mean(R,1));
                else
                    Q = mean(pagemtimes(permute(Psi,[3 2 1]),permute(Psi,[2 3 1])),3);
                    R = mean(pagemtimes(permute(Psi,[3 2 1]),pagemtimes(lambda,permute(Psi,[2 3 1]))),3);
                end
                covPk = ((Q\R)/Q)/sum(Nss-dmax); %Q^-1 R Q^-1
            elseif ny == 1
                W = Nss/sum(Nss);
                Ff = obj.calcKStepPredictor(Hp);
                P = zeros(npar,npar,Hp);
                for j = 0:Hp-1
                    for ns = 1:Nsets
                        Psi = squeeze(de{ns}); 
                        Psi = Psi(~any(isnan(Psi),2),:);
                        X = Psi(1:end-j,:)'*Psi(1+j:end,:)/(length(Psi)-j);
                        P(:,:,j+1) = P(:,:,j+1) + W(ns)*X;
                    end
                    P_H(:,:,j+1) = P(:,:,j+1)*(Ff{1}(1:end-j)*Ff{1}(j+1:end)');
                end
                if Hp > 1 X = 2*sum(P_H(:,:,2:end),3); else X = zeros(npar); end
                P0 = (P_H(:,:,1) + X);
                H_V = P(:,:,1);
                covPk = (H_V\P0)/H_V;
                covPk = lambda*covPk/sum(Nss-dmax);
            end
        end
        
        %% Backcast Ics
        function [e_hat,z0,w_hat,de_hat_dp,z0_dp] = icBackcast(obj,w,ny,dw_dTheta)
            %Computes initial conditions for ResidualFilter CEe = FADEy - FDB f(u,alpha) of a timeseries Model where the input 
            %w is assumed to be the sequence w =  FADEy - FDB f(u,alpha).
            %The overall dimension of the parametervector is assumed to be dim(Theta) = p
            %
            % [e_hat,z0,w_hat,de_hat_dp,z0_dp] = icBackcast(obj,w,ny,dw_dTheta)
            % INPUTS:
            % obj [NonlinearPolyModel]: The Model object
            % w [cell of double vectors or double vector]:  The sequence(s) w. Compute be MA-Filtering w = FADEy - FDB f(u,alpha). If w is a cell
            %                                               the Alg. assumes Nset independet Datasets for which the ICs will be comp.                          	
            % dw_dTheta [opt. cell of N x p matrices]:  Derivatives of sequence w wrt the Parameter vector where N is the number of samles of
            %                                           the corresponding dataset. 
            % ny [pos. int scalar]: Number of Output for which backcasting of IC's will be done
            % OUTPUTS:
            % e_hat [nc x Nset double matrix]: Init. Cond. e(ts-1),...,e(ts-nc) for each of the NSet sequences w, where ts max(na,nb)
            % z0 [nc x Nset double matrix]: Init. Cond. for matlab filter implementation for each of the NSet sequences w. 
            % de_hat_dp [Nset cell of nc x p double matrix]: Init. Cond. de(ts-1),...,e(ts-nc) for each of the NSet sequences w, where ts max(na,nb)
            % z0_dp [Nset cell of nc x p double matrix]: Init. Cond. de_hat_dp for matlab filter implementation.
            % w_hat [nc x Nset double matrix]: Backcasted values of w

            assert(all(obj.Nf(:)==0));
            if nargin <= 2 || isempty(ny); ny = 1; end
            cellin = 1; 
         	nc = obj.Nc(ny);
            ne = obj.Ne(ny);
            nce = nc + ne;
          	ncew = nce + length(obj.PreFilter(ny).den) - 1;
            if ~iscell(w); w = {w}; cellin = 0; end
            if nargin>3 && ~iscell(dw_dTheta); dw_dTheta = {dw_dTheta}; end
            Nset = length(w);
            if obj.C(ny).factorized; [C_,dC] = idModels.util.defactorizePoly(obj.C(ny).val); else; C_ = obj.C(ny).val; end
            if obj.E(ny).factorized; [E_,dE] = idModels.util.defactorizePoly(obj.E(ny).val); else; E_ = obj.E(ny).val; end
            CE = conv(C_,E_);
            CEW = conv(obj.PreFilter(ny).den,CE);
            if nargin == 4; npar = size(dw_dTheta{1},2); end

            % Comp necessary matrices
            Hc = hankel(CEW(2:end));
            Hh = hankel(CEW(end-1:-1:1));
            Hh = Hh(:,end:-1:1);
      
            % Comp d(CEW)
            dCEW_dc = NaN(nc,ncew+1); % C Parameter zeilenweise
            dCEW_de = NaN(ne,ncew+1); % E Parameter zeilenweise
            for k = 1:nc
                if obj.C(ny).factorized
                    dCEW_dc(k,:) = conv(conv(E_,obj.PreFilter(ny).den),dC(k+1,:));
                else
                    dCEW_dc(k,:) = conv(conv(E_,obj.PreFilter(ny).den),[zeros(1,k) 1 zeros(1,nc-k)]);
                end
            end
            dCEW_dc = dCEW_dc(obj.C(ny).free(2:end),:);
            for k = 1:ne
                if obj.E(ny).factorized
                    dCEW_de(k,:) = conv(conv(C_,obj.PreFilter(ny).den),dE(k+1,:));
                else
                    dCEW_de(k,:) = conv(conv(C_,obj.PreFilter(ny).den),[zeros(1,k) 1 zeros(1,ne-k)]);
                end
            end
            dCEW_de = dCEW_de(obj.E(ny).free(2:end),:);

            w_hat = NaN(ncew,Nset); e_hat = w_hat; z0 = w_hat;
            de_hat_dp = cell(Nset,1); z0_dp = de_hat_dp;
            for k = 1:Nset
                % 1. Step: Backwardfiltering of sequence w = CEWe -> e = 1/CEW*w
                [eb_hat, ~] = filter(1,CEW,w{k}(end:-1:1)); %[eb_hat(N) eb_hat(N-1)... eb_hat(ts)]
                %eb_hat = [eb_hat; zeros(ncew-length(eb_hat),1)]; % zeropad if w is too short

                % 2.Step: Backforecast w(k) for k<ts 
                w_hat(:,k) = Hc*eb_hat(end:-1:end-ncew+1); %w_hat(ts-1) ... w_hat(ts-nc)

                % 3. Step: Backforecast e(k) for k<ts 
                e_hat(:,k) = Hh\w_hat(:,k); % e_hat(ts-1) ... e_hat(ts-nc)
                if nargout > 1
                    z0(:,k) = filtic(1,CEW,e_hat(:,k));
                end

                % Comp. Gradient
                if nargin > 3 && nargout > 3
                    % 1. Step: Backwardfiltering of sequence 
                    % dw = d(CEW) = CEW*de + [0 ..., 0 ... 0, e(t+1) ... e(t+nc)] -> de = 1/CEW*(dw - [0 ... 0, 0 ... 0, e(t+1) ... e(t+nc)])
                    X = dw_dTheta{k}(end:-1:1,:); % Reverse sequence
                    ixp = sum([obj.A(ny,:).free obj.B(ny,:).free]) + (1:sum(obj.C(ny).free(2:end)));
                    for np = 1:length(ixp)
                        X(:,ixp(np)) = X(:,ixp(np)) - filter(dCEW_dc(np,:),1,eb_hat);
                    end
                    R2 = zeros(ncew,npar); % needed for step 2        
                    He2 = hankel(eb_hat(end-ncew+1:end)); He2 = He2(:,end:-1:1);
                    R2(:,ixp) = He2*dCEW_dc(:,2:end)';
                    R3 = R2; % needed for step 3
                    He3 = hankel([e_hat(2:end,k)' 0])';
                    R3(:,ixp) = He3*dCEW_dc(:,2:end)';

                    ixp = sum([obj.A(ny,:).free obj.B(ny,:).free obj.C(ny,:).free obj.D(ny,:).free]) + (1:sum(obj.E(ny).free(2:end)));
                    for np = 1:length(ixp)
                        X(:,ixp(np)) = X(:,ixp(np)) - filter(dCEW_de(np,:),1,eb_hat);
                    end
                    R2(:,ixp) = He2*dCEW_de(:,2:end)';
                    R3(:,ixp) = He3*dCEW_de(:,2:end)';
                    deb_hat_dp = filter(1,CEW,X); %[deb(N) deb(N-1)... deb(ts)]

                    % 2. Step Backforecast dw(k) for k<ts
                    dw_hat_dp = Hc*deb_hat_dp(end:-1:end-ncew+1,:) + R2;

                    % 3. Step: Backforecast de(k) for k<ts 
                    de_hat_dp{k} = Hh\(dw_hat_dp - R3);
                    if nargout > 4
                        for l = 1:size(X,2); z0_dp{k}(:,l) = filtic(1,CEW,de_hat_dp{k}(:,l)); end
                    end
                end
            end
            if ~cellin
                if nargin > 3 && nargout > 3
                    de_hat_dp = de_hat_dp{1};
                end
                if nargout > 4
                    z0_dp = z0_dp{1};
                end
            end
        end
        
        %% Residual computation
        function [e, de] = ts_resid(obj,y,u,opt)
        % Computes residuals e and their derivatives de/dP wrt the Parametervecor P.
        % See getPvec and setPvec functions for definition of Parametervector.
        %
        % [e, de] = ts_resid(obj,y,u,opt)
        % y [(cell of) N x 1 double]: Columnvectors of output measurements. Each cell corresponts to a separate dataset.
        % u [(cell of) N x nu double]: Columnvectors of input measurements. Each cell corresponts to a separate dataset.
        % y [(cell of) N x 1 double]: Columnvectors of output measurements. Each cell corresponts to a separate dataset.
        % opt [opt. string]: options structure

        assert(all(obj.Nf(:)==0),'Only Systems with Nf = 0 supported yet!');
        
        %% Init
        cellin = 1; if ~iscell(y); y = {y}; cellin = 0; if ~iscell(u); u = {u}; end; end
        Ns = length(y);
        Nss = cellfun(@length,y);
        N = cellfun(@length,y); 
        nb = obj.Nb; % Order of free part of B
        nc = obj.Nc; nd = obj.Nd;
        nv = length(nb);
        ny = obj.OutputDimension;
        npar = length(obj.getPvec());
        n_alpha = sum(arrayfun(@(fi) sum(fi.free),obj.InputNonlinearity)); %number of free parameters of InputNonlinearity
        dmax = max(obj.getDmax);

        %% Evaluate InputNonlinearity and find indexvector of alpha wrt Pvec
        if strcmpi(opt.Feedback,'sim') && any(obj.HasOutputFeedback)
            assert(obj.OutputDimension == 1);
            obj.UpdateFlag = true;
            obj.updateSs;
            obj.UpdateFlag = false;
            yp = obj.calcPredictions(y,u,'Hp',opt.Hp,'Nmax','auto','UseGls',false,'Warnings',false); %calc residuals init observer with standard LS
            yp = arrayfun(@(ns) [y{ns}(1:opt.Hp-1); yp{ns}(~isnan(yp{ns}(:,:,end)),:,end)],(1:Ns)','UniformOutput',false);
            [fu, dfu] = obj.evalInputNonlinearity(u,yp); % eval inputNonlin
        else  
            [fu, dfu] = obj.evalInputNonlinearity(u,y); % eval inputNonlin
        end
        [n_alpha_i,free_i,ix_alpha_free] = locGetInputNonlinearityIndices();
        
        %% REMOVE OUTPUT OFFSET
        if ~all(obj.OutputOffset==0)
            y = cellfun(@(yi) yi - obj.OutputOffset',y,'UniformOutput',0);
        end

        %% Malloc
        e = arrayfun(@(nsi) NaN(nsi,ny),Nss,'UniformOutput',false); 
        de = arrayfun(@(nsi) zeros(nsi,ny,npar),Nss,'UniformOutput',false);

        %% Calc residuals for each Output
        ixs = 1; % counter for parameters
        for no = 1:ny 
            %% Calculate some Polynomials
            A = arrayfun(@(ai) obj.getCoeff(ai),obj.A(no,:),'UniformOutput',false);
            B = arrayfun(@(bi) obj.getCoeff(bi),obj.B(no,:),'UniformOutput',false);
            C = obj.getCoeff(obj.C(no));
            D = obj.getCoeff(obj.D(no));
            E = obj.getCoeff(obj.E(no));
            L = conv(obj.PreFilter(no).den,conv(C,E)); nl = length(L)-1;      
            AZp = cellfun(@(ai) conv(ai,obj.PreFilter(no).num),A,'UniformOutput',false);
            BZp = cellfun(@(bi) conv(bi,obj.PreFilter(no).num),B,'UniformOutput',false);
            [~,ix_free,~] = obj.getPvec(no);
            npar_no = sum(ix_free);
            ixe = ixs + npar_no - 1 - n_alpha;

            %% Handle ic
            switch lower(opt.Ic)
                case 'zero'
                    ic_meth = 0;
                case 'backcast'
                    ic_meth = 1;
                case 'ls' 
                    ic_meth = 2;
                    if nl>0
                        if isequal(opt.LsInitSamples,Inf)
                            Ninit = max(N)-dmax;
                        elseif strcmpi(opt.LsInitSamples,'Hp')
                            Ninit = max(opt.Hp,nl);
                        elseif strcmpi(opt.LsInitSamples,'auto')
                            Ninit = 10*nl;
                        elseif opt.LsInitSamples>0 && mod(opt.LsInitSamples,1) == 0
                            Ninit = opt.LsInitSamples;
                        end
                        PSI = idModels.alg.lti.getPsi(L,max(N)-dmax); %we need the full PSI because we use it for simulation e = e + PSI*e0 
                    end
                otherwise
                    error('Wrong value of opt.Ic!');
            end

            for ns = 1:Ns
                %% Comp. w: CENp e = w = ADEZp y - BDZp/F*g(u,y) !!
                AZpy_sum = sum(cell2mat(arrayfun(@(no2) filter(AZp{no2},1,y{ns}(:,no2)),1:ny,'UniformOutput',false)),2);
                w = filter(conv(E,D),1,AZpy_sum);
                if nv > 0
                    BZpu_sum = sum(cell2mat(arrayfun(@(ix_u) filter(BZp{ix_u},1,fu{ns}(:,ix_u)),1:nv,'UniformOutput',0)),2);
                    w = w - filter(D,1,BZpu_sum);
                end

                %% Comp e
                if nl==0
                    e{ns}(dmax+1:end,no) = w(dmax+1:end);
                else
                    if ic_meth < 2
                        if ic_meth == 1 % Backcast ic
                            [e_hat,z0] = obj.icBackcast(w(dmax+1:end),no);
                        elseif ic_meth == 0
                            e_hat = zeros(nl,1);
                            z0 = filtic(1,L,e_hat)';            
                        end
                        ef = filter(1,L,w(dmax+1:end),z0);
                    elseif ic_meth == 2
                        Ni = min(Ninit,N(ns)-dmax);
                        R = filter(1,L,w(dmax+1:end));
                        if ~opt.LsInitWrn; warning off; end
                        e_hat = PSI(1:Ni,:)\(-R(1:Ni));      
                        if ~opt.LsInitWrn; warning on; end
                        ef = PSI(1:N(ns)-dmax,:)*e_hat + R;
                    end
                    ix = max(1,dmax-nl+1);
                    e{ns}(ix:end,no) = [e_hat(end-min(dmax,nl)+1:end); ef];
                end

                %% Compute. Gradient de_dp
                if nargout > 1
                    % calc dw_dp
                    dw_dp = calc_dw_dp(); 
                    de_dp = NaN(size(dw_dp));
                    if nl==0 % CE = 1 -> No filtering necessary
                        de_dp = dw_dp;
                    else % get Initial conditions for gradient filter
                        if ic_meth < 2
                            if ic_meth == 1 % Backcast ic
                                [~,~,~,de_hat,z0_dp] = obj.icBackcast(w(dmax+1:end),no,dw_dp(dmax+1:end,:));
                            elseif ic_meth == 0 % zero
                                de_hat = zeros(nl,sum(ix_free));
                                for l = 1:size(de_hat,2); z0_dp(:,l) = filtic(1,L,de_hat(:,l)); end
                            end
                            def = filter(1,L,dw_dp(dmax+1:end,:),z0_dp);
                        elseif ic_meth == 2
                            Ni = min(Ninit,N(ns)-dmax);
                            dR = filter(1,L,dw_dp(dmax+1:end,:)); % forced part (ic = 0)
                            if ~opt.LsInitWrn; warning off; end
                            de_hat = PSI(1:Ni,:)\(-dR(1:Ni,:));
                            if ~opt.LsInitWrn; warning on; end
                            def = PSI(1:N(ns)-dmax,:)*de_hat + dR;
                        end
                        de_dp(ix:end,:) = [de_hat(end-min(dmax,nl)+1:end,:); def];
                    end
                    de{ns}(:,no,ixs:ixe) = de_dp(:,1:end-n_alpha); % LTI-System
                    de{ns}(:,no,end-n_alpha+1:end) = de_dp(:,end-n_alpha+1:end); %InputNonlin Gradients
                end
            end
            ixs = ixe + 1;
        end
        
        %% Set NoiseVariance
        if strcmpi(opt.Feedback,'sim') && any(obj.HasOutputFeedback)
            obj.NoiseVariance = obj.calcCovariance(e,npar);
        end

        %% calc k-(multistep)-residuals
        if any(opt.Hp > 1)
            if nargout > 1
                [e,de] = ts_resid_multistep(obj,e,de,opt);
            else
                e = ts_resid_multistep(obj,e,de,opt);
            end
        end

        %% set NaNs for unavailable samples
        for no = 1:ny
            if opt.MultiStep == 0 || all(opt.Hp(no) == 1)
                Hp = opt.Hp(no);
            else
                Hp = 1:opt.Hp(no);
            end
            for ns = 1:Ns
                ix = 1;
                for hp = Hp
                    e{ns}(1:min(dmax+hp-1,N(ns)),no,ix) = NaN;
                    if nargout > 1
                        de{ns}(1:min(dmax+hp-1,N(ns)),no,:,ix) = NaN;
                    end
                    ix = ix + 1;
                end
            end
        end
        if ~cellin; e = e{1}; if nargout>1; de = de{1}; end; end

        %% Inline Function for Computing dw_dp
        function  dw_dp = calc_dw_dp()
            % comp. dw_dp such that: L*de_dp  = dw_dp (parameters columnwise)
            dw_dp = NaN(N(ns),npar_no); %[ai bi ci di ei alphai]
            ixp = 1; % Laufindex für Parameter
            % A
            yf = filter(conv(D,conv(E,obj.PreFilter(no).num)),1,y{ns}); 
            for no2 = 1:ny
                if obj.A(no,no2).factorized; [~,dA] = idModels.util.defactorizePoly(obj.A(no,no2).val); end
                for np = find(obj.A(no,no2).free)
                    if obj.A(no,no2).factorized 
                        dw_dp(:,ixp) = filter(dA(np,:),1,yf(:,no2));
                    else 
                        dw_dp(np:end,ixp) = yf(1:end-np+1,no2);
                    end
                    ixp = ixp + 1;    
                end
            end
            % B
            fuf = -filter(conv(obj.PreFilter(no).num,D),1,fu{ns}); 
            for ni = 1:nv
                if obj.B(no,ni).factorized; [~,dB] = idModels.util.defactorizePoly(obj.B(no,ni).val); end
                for np = find(obj.B(no,ni).free)-1 
                    if obj.B(no,ni).factorized
                        dw_dp(:,ixp) = filter(dB(np+1,:),1,fuf(:,ni));
                    else
                        dw_dp(np+1:end,ixp) = fuf(1:end-np,ni); 
                    end
                    ixp = ixp + 1; 
                end
            end
            % C
            if nc(no)>0
                if obj.C(no).factorized; [~,dC] = idModels.util.defactorizePoly(obj.C(no).val); end
                nbc = max(nl-dmax,0); % Number of backcasted samples prior t = 1
                ee = [e_hat(1:nbc); e{ns}(:,no)]; % Add additional  backcasted samples
                efc = filter(conv(E,obj.PreFilter(no).den),1,ee);
                for np = find(obj.C(no).free) 
                    if obj.C(no).factorized
                        efci = filter(dC(np,:),1,efc); 
                        dw_dp(:,ixp) = - efci(nbc+1:end); 
                    else
                        dw_dp(max(np-1-nbc,0)+1:end,ixp) = - efc(1+max(0,nbc-np+1):end-np+1); 
                    end
                    ixp = ixp + 1;
                end
            end
            % D
            if nd(no)>0
                if obj.D(no).factorized; [~,dD] = idModels.util.defactorizePoly(obj.D(no).val); end
                yf = filter(E,1,AZpy_sum);
                for np = find(obj.D(no).free) 
                    if obj.D(no).factorized
                        dw_dp(:,ixp) = filter(dD(np,:),1,yf);
                        if nv > 0
                            dw_dp(:,ixp) = dw_dp(:,ixp) - filter(dD(np,:),1,BZpu_sum);
                        end
                    else
                        dw_dp(np:end,ixp) = yf(1:end-np+1);
                        if nv > 0
                            dw_dp(np:end,ixp) = dw_dp(np:end,ixp) - BZpu_sum(1:end-np+1);
                        end
                    end
                    ixp = ixp + 1; 
                end
            end
            % E
            if nv>0 && obj.Ne(no)>0
                if obj.E(no).factorized; [~,dE] = idModels.util.defactorizePoly(obj.E(no).val); end
                yf = filter(D,1,AZpy_sum);
                efe = filter(conv(C,obj.PreFilter(no).den),1,e{ns}(:,no));
                for np = find(obj.E(no).free)
                    if obj.E(no).factorized
                        dw_dp(:,ixp) = filter(dE(np,:),1,yf) - filter(dE(np,:),1,efe);
                    else
                        dw_dp(np:end,ixp) = yf(1:end-np+1) - efe(1:end-np+1); 
                    end
                    ixp = ixp + 1; 
                end
            end
            % Alpha
            dw_dp(:,ixp:end) = zeros(N(ns),n_alpha); % init dw/dalpha
            for ni = 1:nv
                if size(dfu{ns,ni},2) == n_alpha_i(ni) % function returns deriv wrt all Parameters -> comp. Gradient for free parameters only!
                    x = dfu{ns,ni}(:,free_i{ni});
                elseif size(dfu{ns,ni},2) == sum(free_i{ni}) % function returns deriv wrt free parameter vector
                    x = dfu{ns,ni};
                else 
                    error('Incompatible Dimensions!');
                end
                    
                if ~isempty(x)
                    b_df_dalpha = -filter(conv(B{1,ni},conv(obj.PreFilter(no).num,D)),1,x); 
                    dw_dp(:,ix_alpha_free{ni}+ixp-1) = dw_dp(:,ix_alpha_free{ni}+ixp-1) + b_df_dalpha;
                end
            end
        end
        
        function [n_alpha_i,free_i,ix_alpha_free] = locGetInputNonlinearityIndices()
            n_alpha_i = NaN(1,nv);
            ix_alpha = cell(1,nv); ix_fp_i = cell(1,nv); free_i = cell(1,nv);
            np_alpha = 1;
            for nf = 1:nv
                [p_i,ix_fp_i{nf},free_i{nf}] = obj.getFiParams(nf);
                n_alpha_i(nf) = length(p_i); % length of full parametervector parameters
                ix_alpha{nf} = NaN(1,n_alpha_i(nf)); % indexes wrt alpha
                for npi = 1:n_alpha_i(nf)
                    if free_i{nf}(npi) && ix_fp_i{nf}(1,npi) == nf  % numeric free parameters
                        ix_alpha{nf}(npi) = np_alpha;
                        np_alpha = np_alpha + 1;
                    end
                end            
            end
            for nf = 1:nv 
                for npi = 1:n_alpha_i(nf)
                    if ix_fp_i{nf}(1,npi) ~= nf % referenced parameter
                        ix_alpha{nf}(npi) = ix_alpha{ix_fp_i{nf}(1,npi)}(ix_fp_i{nf}(2,npi));
                    end
                end
            end
            ix_alpha_free = cellfun(@(x) x(~isnan(x)) ,ix_alpha,'UniformOutput',false);
        end
        end
        
        %% Calculate cnfidence bounds for Bode diagram
        function [confBodeMagG, confBodeMagH, G, H, confPhaseG, confPhaseH] = calcConfBode(obj,w,confidence)
            %Calculate Confidence of G and H.
            %
            % [confBodeMagG, confBodeMagH, G, H] = calcConfBode(obj,w,confidence)
            % covG [ny x nu cell of nw x 1 double]: Confidence bound of bode magnitude of G in dB 
            % covH [nw x ny double]: Confidence bound of bode magnitude of H in dB 
            % G [ny x nu cell of nw x 1 complex double]: Complex frquency response of G 
            % H [nw x ny complex double]: Complex frquency response of H
            % obj [idModels.NsfPolyModel]: Instance of object.
            % w [nw x 1 pos. double]: Frequency range 0...pi
            % confidence [opt. pos scalar double 0...1]: Confidence bound (Def. 0.99)

            %% Init
            covPk = obj.Info.CovP;
            Conf = true;
            if ~isnumeric(covPk) || any(isnan(covPk(:))) || ~isdiag(obj.Na) % Check weather supported
                if ~isdiag(obj.Na)
                    warning('Couldnt compute confidence bounds because systems with non-diagonal A matrix arent supported yet!');
                else
                    warning('Couldnt compute confidence bounds because parametercovariance (obj.Info.CovP) is not available!');
                end
                Conf = false;
            end
            
            w = obj.Ts*w(:);
            if nargin<3; confidence = .99; end
            assert(confidence>0 && confidence <1);
            npar = size(covPk,1);
            assert(all(obj.Nf(:) == 0),'Models with F polynomial unsupported yet!');
            [nOut, nIn] = size(obj.B);
            G = cell(nOut,nIn);
            confBodeMagG = G;
            confPhaseG = G;
            H = cell(nOut,1);
         	confBodeMagH = H;
            confPhaseH = H;

            %% Calc covG
            for no = 1:nOut
                if no == 1
                    ny0 = 0;
                else
                    ny0 = length(obj.getPvec(no-1));
                end
                if obj.A(no,no).factorized; [A_,dA] = idModels.util.defactorizePoly(obj.A(no,no).val); else; A_ = obj.A(no,no).val; end
                if obj.E(no).factorized; [E_,dE] = idModels.util.defactorizePoly(obj.E(no).val); else; E_ = obj.E(no).val; end
                AE = conv(A_,E_);            
                np_CD = sum(obj.C(no).free) + sum(obj.D(no).free);
                AE2 = conv(AE,AE); 
                for ni = 1:nIn
                    np = 1 + ny0;
                    if obj.B(no,ni).factorized; [B_,dB] = idModels.util.defactorizePoly(obj.B(no,ni).val); else; B_ = obj.B(no,ni).val; end
                    G{no,ni} = freqz(B_,AE,w);
                    if Conf
                        dG = zeros(length(w),npar);
                        for na = find(obj.A(no,no).free)-1
                            if obj.A(no,no).factorized
                                num = - conv(B_,conv(dA(na+1,:),E_));
                            else
                                num = - conv(B_,[zeros(1,na) E_]);
                            end
                            dG(:,np) = freqz(num,AE2,w);
                            np = np + 1;
                        end
                        % A is restricted to be diagonal -> otherwise check here!
                        np = np + sum(arrayfun(@(bi) sum(bi.free),obj.B(no,1:ni-1)));
                        for nb = find(obj.B(no,ni).free)-1
                            if obj.B(no,ni).factorized
                                dG(:,np) = freqz(dB(nb+1,:),AE,w);
                            else
                                dG(:,np) = freqz([zeros(1,nb) 1],AE,w);
                            end
                            np = np + 1;
                        end
                        np = np + sum(arrayfun(@(bi) sum(bi.free),obj.B(no,ni+1:end))) + np_CD;
                        for ne = find(obj.E(1).free)-1
                            if obj.E.factorized
                                num = - conv(B_,conv(dE(ne+1,:),A_));
                            else
                                num = - conv(B_,[zeros(1,ne) A_]);
                            end
                            dG(:,np) = freqz(num,AE2,w);
                            np = np + 1;
                        end
                        % Calc covariance/confidence bounds of Bode Magnitude plot
                        df = 10*(dG.*conj(G{no,ni}) + conj(dG).*G{no,ni})./(log(10)*abs(G{no,ni}).^2);
                        confBodeMagG{no,ni} = util.norminv((confidence+1)/2,zeros(length(w),1),sqrt(diag(df*covPk*df')));
                        % Calc covariance/confidence bounds of phase
                        df = (imag(dG).*real(G{no,ni}) - imag(G{no,ni}).*real(dG))./((1+(imag(G{no,ni})./real(G{no,ni})).^2).*real(G{no,ni}).^2);
                        confPhaseG{no,ni} = util.norminv((confidence+1)/2,zeros(length(w),1),sqrt(diag(df*covPk*df')));
                    else
                        confBodeMagG{no,ni} = NaN(length(w),1);
                        confPhaseG{no,ni} = NaN(length(w),1);
                    end
                end
                %% calc covH
                if nargout > 1
                    if obj.C(no).factorized; [C_,dC] = idModels.util.defactorizePoly(obj.C(no).val); else; C_ = obj.C(no).val; end
                    if obj.D(no).factorized; [D_,dD] = idModels.util.defactorizePoly(obj.D(no).val); else; D_ = obj.D(no).val; end
                    AD = conv(A_,D_);
                    H_fixed = freqz(obj.PreFilter(no).den,obj.PreFilter(no).num,w);
                    H{no,1} = freqz(conv(C_,obj.PreFilter(no).den),conv(AD,obj.PreFilter(no).num),w);
                    if Conf
                        np = 1 + ny0; 
                        AD2 = conv(AD,AD); 
                        Hfree = freqz(C_,AD,w);
                        dH = zeros(length(w),npar);
                        for na = find(obj.A(no,no).free)-1
                            if obj.A(no,no).factorized
                                num = - conv(C_,conv(dA(na+1,:),D_));
                            else
                                num = - conv(C_,[zeros(1,na) D_]);
                            end
                            dH(:,np) = freqz(num,AD2,w);
                            np = np + 1;
                        end
                        np = np + sum(arrayfun(@(bi) sum(bi.free),obj.B(no,:)));
                        for nc = find(obj.C(no).free)-1
                            if obj.C(no).factorized
                                dH(:,np) = freqz(dC(nc+1,:),AD,w);
                            else
                                dH(:,np) = freqz([zeros(1,nc) 1],AD,w);
                            end
                            np = np + 1;
                        end
                        for nd = find(obj.D(no).free)-1
                            if obj.D(no).factorized
                                num = - conv(C_,conv(dD(nd+1,:),A_));
                            else
                                num = - conv(C_,[zeros(1,nd) A_]);
                            end
                            dH(:,np) = freqz(num,AD2,w);
                            np = np + 1;
                        end
                        % Calc covariance/confidence bounds of Bode Magnitude
                        df = 10*(dH.*conj(Hfree) + conj(dH).*Hfree)./(log(10)*abs(Hfree).^2);
                        confBodeMagH{no,1} = util.norminv((confidence+1)/2,zeros(length(w),1),diag(df*covPk*df'));
                        % Calc covariance/confidence bounds of phase
                        df = (real(H{no}).*(real(H_fixed).*imag(dH) + imag(H_fixed).*real(dH)) - ...
                              imag(H{no}).*(real(H_fixed).*real(dH) + imag(H_fixed).*imag(dH)))./ ...
                              ((1+(imag(H{no})./real(H{no})).^2).*real(H{no}).^2);
                        confPhaseH{no,1} = util.norminv((confidence+1)/2,zeros(length(w),1),sqrt(diag(df*covPk*df')));
                    else
                        confBodeMagH{no,1} = NaN(length(w),1);
                        confPhaseH{no,1} = NaN(length(w),1);
                    end
                end
            end
        end

    end
    
    %% private hidden
    methods (Access = private, Hidden = true)
        %% Get Max Delay of Model 
        function dmax = getDmax(obj)
            assert(all(obj.Nf(:)==0),'F-Poly is unsupported yet');
            if size(obj.B,2)>0
                delay = arrayfun(@(no) max(obj.Na(no,:)'+obj.Nd+obj.Ne+length(obj.PreFilter(no).num)-1),(1:obj.OutputDimension)');
                delay =  [delay arrayfun(@(no) max(obj.Nb(no,:)+obj.Nk(no,:)+obj.Nd(no)),(1:obj.OutputDimension)')]; 
            else
                delay = arrayfun(@(no) max(obj.Na(no,:)'+obj.Nd+length(obj.PreFilter(no).num)-1),(1:obj.OutputDimension)');
            end
            dmax = max(delay,[],2);
        end
        
        %% Calc Multistep residuals from 1 step residuals
        function [e,de_dp] = ts_resid_multistep(obj,e,de_dp,opt)
            % NOTE: Function assumes diagonal F, which is giuven if A is daigonal
            assert(isdiag(obj.Na));
            Nsets = size(e,1);
            ny = obj.OutputDimension;
          	npar = size(de_dp{1,1},3);
            n_alpha = sum(arrayfun(@(fi) sum(fi.free),obj.InputNonlinearity)); %number of free parameters of InputNonlinearity
            [Fk,~,dF_da,dF_dc,dF_dd] = obj.calcKStepPredictor(opt.Hp);
            ixy = 1;
            for no = 1:ny
                dF_da_no = dF_da{no,no}(obj.A(no,no).free(2:end),:);
                dF_dc_no = dF_dc{no,no}(obj.C(no).free(2:end),:);
                dF_dd_no = dF_dd{no,no}(obj.D(no).free(2:end),:);
                nay =  sum(arrayfun(@(Ai) sum(Ai.free),obj.A(no,1:no)));
                ixa = ixy-1 + (1:nay);
                nby = sum(arrayfun(@(Bi) sum(Bi.free),obj.B(no,:)));
                ncy = sum(obj.C(no).free);
                ixc = ixy-1 + nay + nby + (1:ncy);
                ndy = sum(obj.D(no).free);
                ixd = ixy-1 + nay + nby + ncy + (1:ndy);
                for ns = 1:Nsets
                    %% Treat Hp-step/Hp-Multistep criterion
                    if ~opt.MultiStep %k-Step
                        if nargout > 1
                            % e_Hp = F e_1
                            % de_Hp/dp_i = dF/dp_i e_1 + F de_1/dp_i
                            de_dp{ns}(:,no,:) = filter(Fk{no,no},1,squeeze(de_dp{ns}(:,no,:))); % F de_1/dp_i (Diagonal F)
                            for ixp = 1:sum(obj.A(no,no).free) % Assumes diagonal A -> diagonal F
                                ef = filter(dF_da_no(ixp,:),1,e{ns}(:,no));
                                de_dp{ns}(:,no,ixa(ixp)) = de_dp{ns}(:,no,ixa(ixp)) + ef; 
                            end
                            for ixp = 1:sum(obj.C(no).free)
                                ef = filter(dF_dc_no(ixp,:),1,e{ns}(:,no));
                                de_dp{ns}(:,no,ixc(ixp)) = de_dp{ns}(:,no,ixc(ixp)) + ef; 
                            end
                            for ixp = 1:sum(obj.D(no).free)
                                ef = filter(dF_dd_no(ixp,:),1,e{ns}(:,no));
                                de_dp{ns}(:,no,ixd(ixp)) = de_dp{ns}(:,no,ixd(ixp)) + ef; 
                            end
                        end
                        e{ns}(:,no) = filter(Fk{no,no},1,e{ns}(:,no)); % F de_1/dp_i
                    else %multistep
                        Ns = size(e{ns},1);
                        em = NaN(Ns,opt.Hp(no));  %[e(t|t-1); e(t|t-2); ... e(t|t-k)]
                        em(:,1) = e{ns}(:,no,1); 
                        if nargout > 1
                            dem = NaN(Ns,npar,opt.Hp(no)); % Ns x npar x Hp
                            dem(:,:,1) = de_dp{ns}(:,no,:,1);
                        end
                        for j = 2:opt.Hp(no)
                            %e_k = F[0:k-1]e_1 = e_(k-1) + f_(k-1)e_1(t-k+1)
                            %de_k = de_(k-1) + F_(k-1)de_1 dF_(k-1)e(k-1)  +  (*)
                            em(:,j) = [ NaN(j-1,1); 
                                    	em(j:end,j-1) + Fk{no,no}(j)*e{ns}(1:end-j+1,no,1)]; 
                            if nargout > 1 % compute "uncorrected" gradients iteratively (last two terms of (*))
                                dem(:,:,j) = [  NaN(j-1,npar); 
                                                dem(j:end,:,j-1) + Fk{no,no}(j)*squeeze(de_dp{ns}(1:end-j+1,no,:,1))];
                            end
                        end
                        if nargout > 1 % now correct gradients
                            efa = zeros(Ns,sum(obj.A(no,no).free)); 
                            efc = zeros(Ns,sum(obj.C(no).free));
                            efd = zeros(Ns,sum(obj.D(no).free));
                            for j = 2:opt.Hp(no) 
                                if size(efa,2) > 0   %na(no,no) > 0
                                    efa(j:end,:) = efa(j:end,:) + dF_da_no(:,j)'.*e{ns}(1:end-j+1,no,1);
                                    dem(j:end,ixa,j) = dem(j:end,ixa,j) + efa(j:end,:); 
                                end
                                if size(efc,2) > 0  %nc(no) > 0
                                    efc(j:end,:) = efc(j:end,:) + dF_dc_no(:,j)'.*e{ns}(1:end-j+1,no,1);
                                    dem(j:end,ixc,j) = dem(j:end,ixc,j) + efc(j:end,:); 
                                end
                                if size(efd,2) > 0  %nd(no) > 0
                                    efd(j:end,:) = efd(j:end,:) + dF_dd_no(:,j)'.*e{ns}(1:end-j+1,no,1);
                                    dem(j:end,ixd,j) = dem(j:end,ixd,j) + efd(j:end,:); 
                                end
                            end
                            de_dp{ns}(:,no,:,1:max(opt.Hp)) = NaN;
                            de_dp{ns}(:,no,:,1:opt.Hp(no)) = dem;
                        end
                        e{ns}(:,no,1:max(opt.Hp)) = NaN; 
                        e{ns}(:,no,1:opt.Hp(no)) = em;
                    end
                end
                ixy = ixy + length(obj.getPvec(ny)) -n_alpha;
            end
        end
        
        %% GET PARAMETER VECTOR
        function [freeP,ixFree,allP,A,b,Pnames] = getPvec(obj,ny)
        % Returns the Parameter Vector used for optimization procedure. 
        % The parameters are in the following order: [A_11 ... A_1ny  B_11 ... B_1nu C1 D1 E1 F_11 ... F_1nu
        %                                             ...
        %                                             A_ny1 ... A_nyny B_ny1 ... B_nynu Cny Dny Eny F_ny1 ... F_nynu
        %                                             alpha1 ... alpha_nu]
        %
        % [freeP,ixFree,allP] = getPvec(obj)
        % freeP [double vector]: The free Parametervector Theta.
        % allP [double vector]: The full Parametervector. 
        % ixFree [logical]: freeP = allP(ixFree)
        % ny [opt pos int]: if supplied only the subparametervector of output ny will be returned
        
        allP = []; ixFree = []; Pnames = {};
        if nargin == 1; ny = 1:obj.OutputDimension; end
        for no = ny
            allP = [allP;
                    cell2mat(arrayfun(@(a) a.val',obj.A(no,:)','UniformOutput',0)); ...
                    cell2mat(arrayfun(@(b) b.val',obj.B(no,:)','UniformOutput',0)); ...
                    obj.C(no).val'
                    obj.D(no).val'
                    cell2mat(arrayfun(@(e) e.val',obj.E(no,:)','UniformOutput',0)); ...
                    cell2mat(arrayfun(@(f) f.val',obj.F(no,:)','UniformOutput',0))];
            ixFree = [  ixFree;
                        cell2mat(arrayfun(@(a) double(a.free'),obj.A(no,:)','UniformOutput',0)); ...
                        cell2mat(arrayfun(@(b) double(b.free'),obj.B(no,:)','UniformOutput',0)); ...
                        double(obj.C(no).free')
                        double(obj.D(no).free')
                        cell2mat(arrayfun(@(e) double(e.free'),obj.E(no,:)','UniformOutput',0)); ...
                        cell2mat(arrayfun(@(f) f.free',obj.F(no,:)','UniformOutput',0))];
          	if nargout > 5
                Pnames = [  Pnames;
                            locGetPnames('A',no);
                            locGetPnames('B',no);
                            locGetPnames('C',no);
                            locGetPnames('D',no);
                            locGetPnames('E',no);
                            locGetPnames('F',no)];
            end
        end
        allP = [allP; cell2mat(arrayfun(@(x) [x.parameters{cellfun(@isnumeric,x.parameters)}]',obj.InputNonlinearity,'UniformOutput',0))];
        ixFree = logical([ixFree; double(cell2mat(arrayfun(@(x) x.free',obj.InputNonlinearity,'UniformOutput',0)))]);
        freeP = allP(ixFree);
        if nargout > 5
            for nv = 1:length(obj.InputNonlinearity)
                Eta = []; 
                if ~isempty(obj.InputNonlinearity(nv).parameters)
                    Eta = {obj.InputNonlinearity(nv).parameters{cellfun(@isnumeric,obj.InputNonlinearity(nv).parameters)}}; 
                end
                for np = 1:length(Eta)
                    Pnames = [  Pnames; 
                                ['Eta_' num2str(nv) num2str(np)]];
                end         
            end
            Pnames = Pnames(ixFree);
        end
        
        % return Matrices for constraints
    	if nargout > 3
            A = []; b = [];
            for ni = 1:length(obj.InputNonlinearity)
                Ai = obj.InputNonlinearity(ni).A;
                A = blkdiag(A,Ai);
                b = [b; obj.InputNonlinearity(ni).b];
            end
            np_f = sum(arrayfun(@(fi) length(fi.free),obj.InputNonlinearity)); % number of params of Inputnonlinearities            
            np_lti = length(allP) - np_f;
            Ncon = size(A,1);
            A = [zeros(Ncon,np_lti) A];
            A = A(:,ixFree);
            ix_isConstraint = any(~(A==0),2);
            A = A(ix_isConstraint,:); % remove Constraints that are 0
            b = b(ix_isConstraint);
        end
        
        function Pnames = locGetPnames(P,no)
           Pnames = {};
           Poly = obj.(P);
           for nu = 1:size(Poly,2)
               if Poly(no,nu).factorized
                    Pnames = [Pnames; ['K' lower(P) '_' num2str(no) num2str(nu)]; arrayfun(@(np) [lower(P) '_' num2str(no) num2str(nu) ',' num2str(np-1) '0'],(2:length(Poly(no,nu).val))','UniformOutput',0)];
               else
                    Pnames = [Pnames; arrayfun(@(np) [lower(P) '_' num2str(no) num2str(nu) ',' num2str(np-1)],(1:length(Poly(no,nu).val))','UniformOutput',0)];
               end
           end
        end
        end
        
        %% SET Parameter vector
        function obj = setPvec(obj,P,StdP)
    	% Sets the free Parametervector P to the object.
        %
        % obj = setPvec(obj,P)
        % P [double vector]: The free Parametervector Theta. The Order of the Parameters is defined in getPvec function.
            np = 1;
            nu = size(obj.Nb,2);
            for no = 1:obj.OutputDimension
                for k = 1:obj.OutputDimension
                    npe = np + sum(obj.A(no,k).free)-1;
                    obj.A(no,k).val(obj.A(no,k).free) = P(np:npe); 
                    if nargin > 2
                        obj.A(no,k).std(obj.A(no,k).free) = StdP(np:npe); 
                    end
                    np = npe + 1;
                end
                for k = 1:nu
                    npe = np + sum(obj.B(no,k).free)-1;
                    obj.B(no,k).val(obj.B(no,k).free) = P(np:npe);
                    if nargin > 2
                        obj.B(no,k).std(obj.B(no,k).free) = StdP(np:npe);
                	end
                    np = npe + 1;
                end
                if obj.Nc(no)>0
                    npe = np + sum(obj.C(no).free)-1;
                    obj.C(no).val(obj.C(no).free) = P(np:npe);
                    if nargin > 2
                        obj.C(no).std(obj.C(no).free) = StdP(np:npe);
                	end
                    np = npe + 1;
                end
                if obj.Nd(no)>0
                    npe = np + sum(obj.D(no).free)-1;
                    obj.D(no).val(obj.D(no).free) = P(np:npe);
                    if nargin > 2
                        obj.D(no).std(obj.D(no).free) = StdP(np:npe);
                	end
                    np = npe + 1;
                end
                if nu>0 && obj.Ne(no)>0
                    npe = np + sum(obj.E(no).free)-1;
                    obj.E(no).val(obj.E(no).free) = P(np:npe);
                    if nargin > 2
                        obj.E(no).std(obj.E(no).free) = StdP(np:npe);
                	end
                    np = npe + 1;
                end
             	for k = 1:nu
                    npe = np + sum(obj.F(no,k).free)-1;
                    obj.F(no,k).val(obj.F(no,k).free) = P(np:npe);
                    if nargin > 2
                        obj.F(no,k).std(obj.F(no,k).free) = StdP(np:npe);
                	end
                    np = npe + 1;
             	end
            end
            for k = 1:length(obj.InputNonlinearity)
                ix_numpar = find(cellfun(@isnumeric,obj.InputNonlinearity(k).parameters));
                ix_freenumpar = ix_numpar(obj.InputNonlinearity(k).free);
                for ix = ix_freenumpar
                    obj.InputNonlinearity(k).parameters{ix} = P(np);
                    if nargin > 2
                        obj.InputNonlinearity(k).std(ix) = StdP(np);
                    end
                    np = np + 1;
                end
            end
        end
        
        %% Change Input indexes
        function obj = changeInputIdx(obj,old_ix,new_ix)
            if ~(old_ix == new_ix)
                assert(obj.OutputDimension == 1);
                lab = obj.InputName;
                Imin = obj.InputMin;
                Imax = obj.InputMax;
                lab{old_ix} = obj.InputName{new_ix};
                lab{new_ix} = obj.InputName{old_ix};
              	Imin(old_ix) = obj.InputMin(new_ix);
                Imin(new_ix) = obj.InputMin(old_ix);
             	Imax(old_ix) = obj.InputMax(new_ix);
                Imax(new_ix) = obj.InputMax(old_ix);
                obj.InputName = lab;
                obj.InputMin = Imin;
                obj.InputMax = Imax;
                for ni = 1:length(obj.InputNonlinearity) 
                    ix1 = find(old_ix == obj.InputNonlinearity(ni).input_idx);
                    ix2 = find(new_ix == obj.InputNonlinearity(ni).input_idx);
                    if ~isempty(ix1)
                        obj.InputNonlinearity(ni).input_idx(ix1) = new_ix;
                    end
                    if ~isempty(ix2)
                        obj.InputNonlinearity(ni).input_idx(ix2) = old_ix;
                    end
                end
            end
        end
        
    	%% Stabilization 
        function obj = stabilize(obj,P,fac)
            if nargin == 2
                fac = 1.1;
            end
            for no = 1:obj.OutputDimension
                obj.(P)(no) = idModels.alg.lti.stabilizePoly(obj.(P)(no),fac);
            end
        end
        
       	%% Plot zeros of Polynomial
        function h = plotZeros(~,p,names)
            h = figure('color','w');
            np = length(p);
            for k = 1:np
                P = p{k}([true; sum(abs(diff(p{k})),2)~=0],:);
                Pr = cell2mat(arrayfun(@(j) roots(P(j,:)),1:size(P,1),'UniformOutput',false))';
                l = linspace(0,2*pi,100);
                subplot(1,np,k); plot(cos(l),sin(l),'-k'); xlabel('Real'); ylabel('Imag'); title(['Roots of ' names{k}]);
                for kk = 1:size(Pr,2)
                    absPr = [0; abs(diff(abs(Pr(:,kk))))];
                    col = cumsum(absPr); col = min(.99,abs(col./col(end)-1));
                    hold on; col = kron(col,[1 1 1]);
                    plot(real(Pr(:,kk)),imag(Pr(:,kk)),'LineStyle','-','Color',[.9 .9 .9]);
                    scatter(real(Pr(:,kk)),imag(Pr(:,kk)),48,col); hold on;
                    scatter(real(Pr(end,kk)),imag(Pr(end,kk)),48,'r'); hold on;
                    lim = [get(gca,'XLim')' get(gca,'YLim')']; lim = [min(lim(1,:)) max(lim(2,:))];
                    set(gca,'XLim',lim,'YLim',lim);  
                end
                grid on; axis equal;
            end
        end
        
        %% Calculate Frequency weighting function
        function W = calcFreqWeight(obj,w)
            %Calculate Frequency Weighting Function of Fit.
            %
            % W = calcFreqWeight(obj,w)
            % w [double vector]:
            % W [N x ny complex double]: The weighting function.
            % obj [idModels.NsfPolyModel]: Instance of model object.

            assert(isvector(w) && isreal(w),'w needs to be a real vector!');
            
            ny = obj.OutputDimension;
            W = NaN(length(w),ny);
            if obj.Info.OptionsUsed.MultiStep == true && obj.Info.OptionsUsed.Hp > 1 
                warning('Cant compute Frequency weighting because multistep criterion unsupported yet!');
                return;
            end
            
            w = obj.Ts*w(:);
            [~,~,H_num,H_den] = calcGH(obj,'AbsorbNoiseVariance',false);
            Ff = obj.calcKStepPredictor(obj.Info.OptionsUsed.Hp);
            
            for no = 1:ny
                num_W = conv(Ff{no,no},H_den{no});
                den_W = H_num{no};
                W(:,no) = freqz(num_W,den_W,w);
            end
        end
                
        %% Comp initial Values
        function obj = initPoly(obj,y,u,opt)
            % Initializes Polynomial Model using LS/IV techniques.

            up = obj.evalInputNonlinearity(u,y); % Evaluate Inputnonlinearities
            y = cellfun(@(yi) bsxfun(@minus,yi,obj.OutputOffset'),y,'UniformOutput',0); % Now Remove OutputOffset
            
            % Preinit
            ny = obj.OutputDimension;
            nu = size(obj.B,2);

            % Treat partial Init
            a = obj.A; b = obj.B; c = obj.C; d = obj.D; e = obj.E; 
            
            % Fix parameters that have initial guess, compute coeff of factorized polynomials
            for nr = 1:ny
                for nc = 1:ny; if a(nr,nc).factorized; a(nr,nc) = localInitFacPoly(a(nr,nc)); end; end
                for nc = 1:nu; if b(nr,nc).factorized; b(nr,nc) = localInitFacPoly(b(nr,nc)); end; end
                if c(nr).factorized; c(nr) = localInitFacPoly(c(nr)); end
                if d(nr).factorized; d(nr) = localInitFacPoly(d(nr)); end
                if e(nr).factorized; e(nr) = localInitFacPoly(e(nr)); end
            end

            % Prefilter
            Pf = cell(ny,1);
            for no = 1:ny
                Pf{no} = {obj.PreFilter(no).num obj.PreFilter(no).den};
            end

            % Assign Outputs
            [A0,B0,C0,D0,E0,~] = idModels.util.polyStruct2cell(obj.A,obj.B,obj.C,obj.D,obj.E,obj.F,'SetNanTo0',true);

            % Start Init
            if strcmpi(opt.InitMethod,'plr') || strcmpi(opt.InitMethod,'arx') || strcmpi(opt.InitMethod,'auto') 
                if strcmpi(opt.InitMethod,'arx')
                    opt.InitIter = 1;
                    isEyeA = all(all(arrayfun(@(ai) isequal(ai.val,1) || isequal(ai.val,0),a)));
                    if isEyeA  %use AR part to initilize E (and D) in BJ and OE case
                        oldA0 = A0;
                        for no = 1:ny; a(no,no) = e(no); end
                    end
                end
                
                %% AR(MA)(X) - the AR(X) Case is handled inside plr algorithm
                [A0, B0, C0, D0, E0]  = idModels.alg.ls.plr(y,up,a,b,[],c,d,e,'Init','ls','Plot',0,'Prefilter',Pf,'Iter',opt.InitIter,'Plot',false);
                if strcmpi(opt.InitMethod,'arx') && isEyeA
                    for no = 1:ny
                        E0{no,1} = A0{no,no};
                        % set free D parameters
                        freeD = find(d(no).free(1:min(length(A0{no,no}),length(D0{no}))));
                        D0{no,1}(freeD) = A0{no,no}(freeD);    
                    end
                    A0 = oldA0;
                end
                
                % PLR didnt converge successfully
                if any(any(cellfun(@(x) any(isnan(x)),A0))) || any(any(cellfun(@(x) any(isnan(x)),B0))) || ...
                        any(any(cellfun(@(x) any(isnan(x)),C0))) || any(any(cellfun(@(x) any(isnan(x)),D0))) || any(any(cellfun(@(x) any(isnan(x)),E0)))
                    opt.InitMethod = 'arx';
                    warning('Some Problems occured while initializing Model with PLR-Method - > Trying %s-Method',opt.InitMethod);
                    obj = initPoly(obj,y,u,opt);
                    return;
                end
            else % case iv 
                %% DETERMINISTIC PART
                if nu > 0 % estimate determinitsic part by using IV Method
                    % If BJ -> init B,E
                    % If AR... -> init B,A
                    if all(obj.Na(:) == 0)  % BJ: Set aa = e
                        a0.val = 0; a0.free = false;
                        for no1 = 1:ny
                            for no2 = 1:ny
                                if no1 == no2
                                    aa(no1,no2) = e(no1);
                                else
                                    aa(no1,no2) = a0;
                                end
                            end
                        end
                    else % ARARMAX or full model
                        aa = a;
                    end
                    [AA,B0,~,v] = idModels.alg.ls.oe_iv(y,up,aa,b,[],'Iter',opt.InitIter,'Plot',0,'EstIc',0,'Prefilter',Pf,...
                                        'StabilizationFactor',opt.StabilizationFactor,'Ls_init_samples','auto'); 

                    if any(any(cellfun(@(x) any(isnan(x)),AA))) || any(any(cellfun(@(x) any(isnan(x)),B0)))
                        opt.InitMethod = 'PLR';
                        warning('Some Problems occured while initializing Model with IV-Method - > Trying %s-Method',opt.InitMethod);
                        obj = initPoly(obj,y,u,opt);
                        return;
                    end

                    % Remove NaNs
                    for ns = 1:length(y)
                        v{ns} = v{ns}(all(~isnan(v{ns}),2),:); % Remove NaNs
                    end

                    % assign 
                    if all(obj.Na(:) == 0)  % BJ
                        E0 = AA;
                    else
                        A0 = AA;
                    end        
                else % Autonomous case
                    v = y;
                end
                %% Now estimate Noise Model with PEM
                if strcmpi(opt.InitMethod,'iv') && any(obj.Nd + obj.Nc > 0) || (nu == 0 && any(obj.Na >0))
                    dum = idModels.NsfPolyModel(obj.Na,[],[],obj.Nc,obj.Nd,'C',c,'D',obj.D); %ARMA
                    if nu > 0 && ~ all(obj.Na(:) == 0) % Fix A if A was identified already
                        for ny1 = 1:ny; for ny2 = 1:ny; dum.A(ny1,ny2).val = A0{ny1,ny2}; dum.A(ny1,ny2).free(1:end) = 0; end; end 
                    end
                    if ~dum.IsParameterized
                        opt.MaxIter = 1e2;
                        opt.Hp = ones(ny,1); 
                        opt.InitMethod = 'arx';
                        opt.EstimateOutputOffset = false; 
                        opt.PlotZeros = {}; 
                        opt.Display = 'iter'; 
                        opt.IntegrateNoise = false(ny,1);
                        opt.CheckGradient = 'off';
                        opt.ForcePolesG = '';
                        opt = rmfield(opt,{'u_norm' 'y_norm' 'Unmatched'});
                        dum.identify(v,opt);
                    end
                    A0 = idModels.util.polyStruct2cell(dum.A);
                    C0 = idModels.util.polyStruct2cell(dum.C);
                    D0 =  idModels.util.polyStruct2cell(dum.D);
                    clear dum;
                end
            end
            
            % Init nonfactorized polynomials only (factorized polys are assumed to have initial guesses already) 
            for nr = 1:ny
                for nc = 1:ny
                    if any(isnan(obj.A(nr,nc).val))
                        if ~obj.A(nr,nc).factorized; obj.A(nr,nc).val = A0{nr,nc}; else; obj.A(nr,nc).val = idModels.util.factorizePoly(A0{nr,nc}); end 
                    end
                end
                for nc = 1:nu 
                    if any(isnan(obj.B(nr,nc).val))
                        if ~obj.B(nr,nc).factorized; obj.B(nr,nc).val = B0{nr,nc}; else; obj.B(nr,nc).val = idModels.util.factorizePoly(B0{nr,nc}); end
                    end
                end
                if any(isnan(obj.C(nr).val)) 
                    if ~obj.C(nr).factorized; obj.C(nr).val = C0{nr}; else; obj.C(nr).val = idModels.util.factorizePoly(C0{nr}); end 
                end
                if any(isnan(obj.D(nr).val)) 
                    if ~obj.D(nr).factorized; obj.D(nr).val = D0{nr}; else; obj.D(nr).val = idModels.util.factorizePoly(D0{nr}); end 
                end
                if any(isnan(obj.E(nr).val)) 
                    if ~obj.E(nr).factorized; obj.E(nr).val = E0{nr}; else; obj.E(nr).val = idModels.util.factorizePoly(E0{nr}); end
                end
            end
            assert(obj.IsParameterized,'Some error occured during initialization!');
            
            function p = localInitFacPoly(p) 
                n = length(p.val)-1;
              	nk = isinf(p.val) & ~p.free; % idx of fixed delays
                hasFreeZero = double((n - sum(nk)) > 0);
                if all(~isnan(p.val))  % Poly is initilized -> set free to false and defactorize 
                  	p.free(1:end) = false; 
                    p.val = idModels.util.defactorizePoly(p.val,obj.ImagTol);
                elseif all(isnan(p.val([false ~nk(2:end)]))) % all poles arent initilized (except delays and gain)
                    p.val = [zeros(1,sum(nk)) p.val(1) NaN(hasFreeZero) zeros(1,n-sum(nk)-1)]; %1st coeff 1
                    p.free = [false(1,sum(nk)) p.free(1) true(hasFreeZero) false(1,n-sum(nk)-1)]; % one coeff is free 
                else 
                    error('Partial initialization of factorized polyniomials is unsupported yet!')
                end
                p.factorized = false;
            end
        end
        
        %% Check Polynomial values
        function P = checkPoly(obj,P,poly)
            % Checks Polynomial P for validity.
            % if flag == 1: A Check:    Checks if first coeff of diag polynomials is 1 and first coeff of offdiag polynomials is 0 
            %                           
            % if flag == 2: B Check
            % if flag == 3: C/D/E/F Check: Make sure that leeding coefficient is 1 and not free
            
            if ~obj.DoChecks
                return;
            end

            for k = 1:size(P,1) 
                for l = 1:size(P,2)
                    if isempty(P(k,l).val) && isempty(P(k,l).factorized); P(k,l).factorized = false; end 
                    assert(all(~isnan(P(k,l).val(~P(k,l).free))),'Fixed parameters need to real numeric values (no NaN or Inf)!');
                    assert(isempty(obj.(poly)) || (obj.(poly)(k,l).factorized==P(k,l).factorized),['The property "factorized" is not allowed to change! Use factorize() and defactorize() method to change the parameterization of the Polynomial ' poly '!']);
                  	assert(isempty(fieldnames(rmfield(P(k,l),{'val' 'free' 'factorized' 'std'}))),'Only "val", "free", "std" and "factorized" are valid fields!');
                    assert((all(isreal(P(k,l).val)) || P(k,l).factorized) && (P(k,l).factorized || all(~isinf(P(k,l).val))) && (isvector(P(k,l).val) || isempty(P(k,l).val)) && (isvector(P(k,l).free) || isempty(P(k,l).free)),'Coefficients of polynomials need to be real vectors!');
                    assert(length(P(k,l).val) == length(P(k,l).free) && length(P(k,l).val) == length(P(k,l).std),'The "val", "std" and the free vector need to have equal length!');
                    assert(all(P(k,l).free==0 | P(k,l).free==1),'The "free" vectors need to be logicals');
                    assert(all(P(k,l).std>=0 | isnan(P(k,l).std)),'The "std" vectors need to be positive');
                    P(k,l).free = logical(P(k,l).free);
                    switch poly
                        case 'A' % A 
                            if P(k,l).factorized
                                if k == l
                                    assert(P(k,l).val(1)==1 & P(k,l).free(1)==0,'Gain of polynomial need to be 1 and not free!');
                                else
                                    assert(obj.getPolyDelay(P(k,l))>0,'Off-diagonal polynomials need to have at least one delay!');
                                end
                            else
                                if isempty(P(k,l).val) %make valid coeff
                                    if k==l
                                        P(k,l).val = 1;
                                    else
                                        P(k,l).val = 0;
                                    end
                                    P(k,l).free = false;
                                    P(k,l).std = NaN;
                                end
                                if k == l
                                    assert(P(k,l).val(1)==1 & P(k,l).free(1)==0,'Leeding coefficients of diagonal polynomials need to be 1 and not free!');
                                else
                                    assert(P(k,l).val(1)==0,'Leeding coefficients of off-diagonal polynomials need to be 0!');
                                end
                            end
                        case 'B' 
                            if isempty(P(k,l).val) %make valid coeff
                                P(k,l).val = 0;
                                P(k,l).free = false;
                                P(k,l).std = NaN;
                            end
                        case {'C' 'D' 'E' 'F'}
                            if isempty(P(k,l).val) %make valid coeff
                                P(k,l).val = 1;
                                P(k,l).free = false;
                                P(k,l).std = NaN;
                            end
                            if P(k,l).factorized
                                assert(P(k,l).val(1)==1 & P(k,l).free(1)==0,'Gain of polynomial need to be 1 and not free!');
                            else
                                assert(P(k,l).val(1)==1 && P(k,l).free(1)==0,'Leeding coefficients of diagonal polynomials need to be 1 and not free!');
                            end
                    end
                    ix = abs(imag(P(k,l).val))<=idModels.NsfPolyModel.ImagTol;
                    P(k,l).val(ix) = real(P(k,l).val(ix));
                end
            end
        end
        
    	function pi = getCoeff(obj,P)
            if P.factorized
                pi = idModels.util.defactorizePoly(P.val,obj.ImagTol);
            else
                pi = P.val;
            end
        end
        
        %% Dereferencing of Parametervector
        function [pi,ix_fp,free] = getFiParams(obj,ix)
            pi = NaN(1,length(obj.InputNonlinearity(ix).parameters));
            ix_str = cellfun(@ischar,obj.InputNonlinearity(ix).parameters);
            for np = find(~ix_str)
                pi(np) = obj.InputNonlinearity(ix).parameters{np};
            end
            if nargout > 2
                free = false(1,length(obj.InputNonlinearity(ix).parameters));
                free(~ix_str) = obj.InputNonlinearity(ix).free;
            end
            ix_fp = [   ix*ones(1,length(pi)) 
                        1:length(pi)];
            for np = find(ix_str)
                str = obj.InputNonlinearity(ix).parameters{np};
                ix_ = strfind(str,'_');
                ix_fp(1,np) = str2double(str(1:ix_-1));
                ix_fp(2,np) = str2double(str(ix_+1:end));
                pi(np) = obj.InputNonlinearity(ix_fp(1,np)).parameters{ix_fp(2,np)};
                if nargout > 2
                    % find corresponding free
                    nrefpar = sum(cellfun(@ischar,obj.InputNonlinearity(ix_fp(1,np)).parameters(1:ix_fp(2,np)-1)));
                    free(np) = obj.InputNonlinearity(ix_fp(1,np)).free(ix_fp(2,np)-nrefpar);
                end
            end
        end
        
    end
    
    %% PRIVATE STATIC 
    methods (Access = private, Hidden = true, Static = true)                 
        function [Pc,Rem] = commonZeros(P,PolyOut)
            if nargin == 1
                PolyOut = false;
            end
            Pc = [];
            AllZeros = cellfun(@roots,P,'UniformOutput',0);
            Rem = AllZeros;
            
            for k = 1:length(AllZeros{1})
                isc = cellfun(@(xi) any(AllZeros{1}(k) == xi),Rem);
                if all(isc)
                    Pc = [Pc AllZeros{1}(k)];
                    Rem = cellfun(@(xi) xi(xi~=AllZeros{1}(k)),Rem,'UniformOutput',0);
                end
            end
            if PolyOut
                Rem = cellfun(@(xi) poly(xi),Rem,'UniformOutput',0);
                Pc = poly(Pc);
            end
        end
        
        function Pstr = setPolyVals(Pstr,P)
           	if isvector(P) && isnumeric(P)
                P = {P};
            end
            if isvector(Pstr) && ~isvector(P) && isdiag(cellfun(@length,P)-1)
                P = arrayfun(@(i) P{i,i},(1:length(P))','UniformOutput',false);
            end
            [nr, nc] = size(Pstr);
            assert(all([nr nc] == size(P)),'Incompatible Dimension for setting poly values!'); 
            for r = 1:nr
                for c = 1:nc
                    assert(isvector(P{r,c}),'The values to be set need to be double vectors!');
                    assert(length(Pstr(r,c).val) == length(P{r,c}),'Incompatible Dimensions!');
                    Pstr(r,c).val = P{r,c};
                    Pstr(r,c).std = NaN(1,length(P{r,c}));
                end 
            end
        end
        
        function n = getPolyDegree(P)
            n = NaN(size(P));
            for k = 1:size(P,1) 
                for l = 1:size(P,2)
                    n(k,l) = length(P(k,l).val) - 1;
                end
            end
        end
        
        function nk = getPolyDelay(P)
        	nk = NaN(size(P));
            for k = 1:size(P,1) 
                for l = 1:size(P,2)
                    if P(k,l).factorized
                        if P(k,l).val(1) == 0 && P(k,l).free(1) == false % Fixed 0 gain 
                            nk(k,l) = length(P(k,l).val); 
                        else
                            nk(k,l) = sum(P(k,l).val([false ~P(k,l).free(2:end)])==Inf); % Number of fixed Zeros at Infinity
                        end
                    else
                        d = find((P(k,l).free==false & P(k,l).val==0)==0,1,'first')-1; % Number of fixed coeff that are 0
                        if ~isempty(d)
                            nk(k,l) = d;
                        else
                            nk(k,l) = length(P(k,l).val);
                        end
                    end
                end
            end 
        end
        
        function h = polyHasZero(P,val,tol)
            if nargin < 3; tol = 1e-8; end
            if ~P.factorized && all(~isnan(P.val))
                rD = roots(P.val);
            elseif P.factorized
                rD = P.val(2:end); 
            else
                rD = [];
            end
            if ~isempty(rD) && any(abs(rD-val) <= tol) 
                h = true;
            else
                h = false;
            end 
        end
    end
    
	%% Static Methdos
    methods (Static = true)
        [Cost,M] = estModel(y,u,varargin)
        
        function obj = loadobj(str)
            if isfield(str,'Ss'); Ssobj = str.Ss ; str = rmfield(str,'Ss'); end
            if isfield(str,'F_nargout'); str = rmfield(str,'F_nargout'); end
            if isfield(str,'DoChecks'); str = rmfield(str,'DoChecks'); end
            if isfield(str,'UpdateFlag'); str = rmfield(str,'UpdateFlag'); end
            obj = idModels.NsfPolyModel(zeros(size(str.A)),zeros(size(str.B)),zeros(size(str.B)),zeros(size(str.C)),zeros(size(str.D)),zeros(size(str.E)),zeros(size(str.F)),...
                                'InputName',str.InputName,'OutputName',str.OutputName,'InputNonlinearity',str.InputNonlinearity);
            obj.DoChecks = false;
            obj = obj.set(str);
            obj.Ss = Ssobj;
            obj.DoChecks = true;
            obj.UpdateFlag = true;
        end
    end
end


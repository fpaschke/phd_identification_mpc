classdef PolynomialModel < idModels.StaticModel
    %POLYNOMIAL Defines a polynomial Model m inputs with orders n1 ... nm
    %   y(u,p) = p0 + 
    %   p11*u1 + p12*u1^2 + ... + p1n1*u1^n1 + ...
    %   ...
    %   pm1*um + pm2*um^2 + ... + p1nm*um^nm. 
    %
    % obj = PolynmomialModel(degree,varargin)
    % obj [idModels.PolynomialModel]: The object to be created
    % varargin: Name Value Pairs specifing the object. You can use any valid property (which is not dependent)
    % of PolynomialModel (or its parent eg IDModel) and pass it to the constructor.
    %   
    % EXAMPLES: 
    %
	% M = idModels.PolynomialModel(3,'InputName',{'u1'},'OutputName',{'y1'},'OutputOffset',10);
    % N = 1e3; U = 10*(rand(N,1)-.5); Y = 20 + 2*U(:,1) + 1*U(:,1).^2 + 1*U(:,1).^3 + 10*randn(N,1); 
    % M.Coefficients.free(1) = false; M.Coefficients.val(1) = 2; %fix linear and quad part
    % M.Coefficients.free(2) = false; M.Coefficients.val(2) = 1; %fix linear and quad part
	% iopt = M.identifyOptions('EstimateOutputOffset',true); 
    % M.identify(Y,U,iopt);
    % M.show('Autocovariance',Y,U,'Significance',.01)
    %
    % M = idModels.PolynomialModel(6);
    % N = 5; rng(0); U = linspace(0,2,N)'; Y = 2 + 2*U + 1*U.^2 + 3*U.^3 + 5*randn(N,1); 
    % M.identify(Y,U,'EstimateOutputOffset',1,'Nconstraints',10,'FirstDerivative',1,'SecondDerivative',1)
    % M.showModel('Data',{Y U})
    %
    % N = 100; rng(0); U = 2*rand(N,3); Y = 2 + 1*U(:,1) + 2*U(:,2) + 3*U(:,3) + 5*randn(N,1); 
    % M1 = idModels.PolynomialModel([2 2 1],'OutputName','PolyModel_1'); M1.identify(Y,U,'Nconstraints',100,'FirstDerivative',0,'SecondDerivative',0,'EstimateOutputOffset',1);
    % M2 = idModels.PolynomialModel([2 2 1],'OutputName','PolyModel_2'); M2.identify(Y,U,'Nconstraints',100,'FirstDerivative',1,'SecondDerivative',0,'EstimateOutputOffset',1);
  	% M3 = idModels.PolynomialModel([2 2 1],'OutputName','PolyModel_3'); M3.identify(Y,U,'Nconstraints',100,'FirstDerivative',1,'SecondDerivative',1,'EstimateOutputOffset',1);
    % M = merge(M1,M2,M3);
    % M.showModel()
    
    properties
        Coefficients            % [ny x nu struct]: struct with fields val [1 x n double] and free [1 x n double] 
    end
    
    properties (Dependent)
    	Degree                  % [ny x nu pos int]: indicating degrees n of corresponding polyniomials.
    end
    
    properties (Dependent, Hidden)
        HasEqualCoefficients    % [ny x nu pos int]: indicating weather coefficients of inputs are equal. [1 2 1; 1 1 1] indicates that input 1 and 3 of output 1 and input 1, 2 and 3 of input 2 have equal coefficients (Used to reduce computational demand for simulation)
    end
    
    methods %public Methods
        function obj = PolynomialModel(degree,varargin)           
            % Make sure to initialize input and outputlabels to valid values first
            [ny,nu] = size(degree);
            inname = arrayfun(@(i) ['u' num2str(i)],1:nu,'UniformOutput',0)'; 
            outname = arrayfun(@(i) ['y' num2str(i)],1:ny,'UniformOutput',0)';
            
            obj = obj@idModels.StaticModel('InputName',inname,'OutputName',outname);         
            obj.Coefficients = obj.initCoefficientsFromDegree(degree);
            obj = obj.set(varargin{:}); %Now set user supplied opt. name value pairs!
        end
        
        %% Prop getters and setters
        % Coefficients
        function set.Coefficients(obj,coeff)
            % Check if valid and set
            isnumeric_flag = arrayfun(@(x) isnumeric(x.val),coeff);
            assert(all(isnumeric_flag(:)),'All Coefficients need to be numeric values');
            assert(all(all(arrayfun(@(c) all(c.free == 0 | c.free == 1) && length(c.free) == length(c.val),coeff))),'"free" needs to be a logical vector with same length ans val!');
            assert(isequal(size(coeff),[obj.OutputDimension obj.InputDimension]),'The size(Coefficients) needs to be [OutputDimension InputDimension]');  
            for ny = 1:obj.OutputDimension 
                for nu = 1:obj.InputDimension
                    coeff(ny,nu).val = coeff(ny,nu).val(:)';
                    coeff(ny,nu).free = logical(coeff(ny,nu).free(:)'); 
                end
            end
            obj.Coefficients = coeff;
        end
        
        function coeff = get.Coefficients(obj)
           coeff = obj.Coefficients; 
        end
        
        %% Dependent Prop getters and setters
        % HasEqualCoefficients
        function ise = get.HasEqualCoefficients(obj)
            % Checks weather the coefficients of each input are equal and returns
            % a array of integers indicating inputs with equal coefficients
            ise = zeros(obj.OutputDimension,obj.InputDimension);
            for no = 1:obj.OutputDimension
                k = 1;
                while any(ise(no,:) == 0)
                    fn = find(ise(no,:) == 0,1,'first');
                    ix = cellfun(@(x) isequaln(x,obj.Coefficients(no,fn).val),{obj.Coefficients(no,:).val});
                    ise(no,ix) = k;
                    k = k + 1;
                end
            end
        end
        
        % Degree
      	function deg = get.Degree(obj)
            deg = arrayfun(@(x) length(x.val),obj.Coefficients);
        end
    
        %% Abstract implementations
        e = identifyModel(obj,y,u,varargin)
      	opt = identifyModelOptions(obj,varargin)
        y = simulateModel(obj,u)
	end %public Methods
    
    %% PRIVATE METHS
    methods (Access = private, Hidden = true)
        function [PHI,dPHI,ddPHI] = getRegressorMatrices(obj,u)
            %GETREGERSSORMATRICES Caluculates Regressorematrices of Polynomial Model with ny outputs and nu inputs.
            %
            % PHI = getRegressorMatrices(obj,u)
            % u [N x nu double Matrix]: Input Data
            % PHI [ny x nu cell of matrices]: Regressormatrix PHI{ny,nu} [ui.^1 ui.^2 ... ui.^deg(i)]
            % dPHI [ny x nu cell of matrices]: 1st derivative of Regressormatrix dPHI{ny,nu}    [ui.^0 2*ui.^1  3*ui.^2 ... deg(i)*ui.^(deg(i)-1)]
            % ddPHI [ny x nu cell of matrices]: 2nd derivative of Regressormatrix ddPHI{ny,nu}  [0     2*ui.^0  6*ui.   ... deg(i)*ui.^(deg(i)-1)]
            
            deg = obj.Degree;
            Ns = size(u,1);
            assert(size(deg,2) == size(u,2),'Length of deg needs to be size(u,2)!');

            % build Regressor Matrix 
            [ny, nu] = size(deg);
            PHI = cell(ny,nu); dPHI = PHI; ddPHI = PHI;
            for no = 1:ny
                for ni = 1:nu
                    PHI{no,ni} = bsxfun(@power,u(:,ni),1:deg(no,ni));
                    if nargout >= 2
                        dPHI{no,ni} = cell2mat(arrayfun(@(k) k*u(:,ni).^(k-1),1:deg(no,ni),'UniformOutput',0));
                    end
                    if nargout >= 3
                        ddPHI{no,ni} = cell2mat(arrayfun(@(k) k*(k-1)*u(:,ni).^(k-2),2:deg(no,ni),'UniformOutput',0));
                        if deg(no,ni) >= 1 ddPHI{no,ni} = [zeros(Ns,1) ddPHI{no,ni}]; end
                    end
                end
            end
        end
        
        function [P,ixfree,allP] = getPvec(obj,un,yn)
      		%GETPVEC Returns Parametervector. 
            if nargin <= 1 un = ones(1,obj.InputDimension); end
            if nargin <= 2 yn = ones(1,obj.OutputDimension); end
            
          	allP = []; ixfree = logical([]);
            for ny = 1:obj.OutputDimension
                for nu = 1:obj.InputDimension
                    allP = [allP; (obj.Coefficients(ny,nu).val.*bsxfun(@power,un(nu),1:obj.Degree(ny,nu)))'/yn(ny)];
                    ixfree = [ixfree; obj.Coefficients(ny,nu).free'];
                end
            end
            P = allP(ixfree);
        end
        
        function obj = setPvec(obj,P,un,yn)
      		%SETPVEC Sets Parametervector. 
            
            if nargin == 1 
               un = ones(1,obj.InputDimension);
               yn = ones(1,obj.OutputDimension);
            end
            
            ixp = 1;
            for ny = 1:obj.OutputDimension
                for nu = 1:obj.InputDimension
                    ixe = ixp + sum(obj.Coefficients(ny,nu).free)-1;
                    obj.Coefficients(ny,nu).val(obj.Coefficients(ny,nu).free) = ...
                        yn(ny)*bsxfun(@power,1/un(nu),find(obj.Coefficients(ny,nu).free))'.*P(ixp:ixe); 
                    ixp = ixe + 1;
                end
            end
        end
    end % private methods
    
    %% private static methods
    methods (Static = true, Access = private, Hidden = true)      
        function coeff = initCoefficientsFromDegree(deg)
        	for ny = 1:size(deg,1)
                for nu = 1:size(deg,2)
                    coeff(ny,nu).val = NaN(deg(ny,nu),1);
                    coeff(ny,nu).free = true(length(coeff(ny,nu).val),1); 
                end
            end
        end
    end % private static methods
	
    %% Static Methdos
    methods (Static = true)
        function obj = loadobj(str)
            deg = cellfun(@length,str.Coefficients);
            obj = idModels.PolynomialModel(deg);
            obj = obj.set(str);
        end
    end
end


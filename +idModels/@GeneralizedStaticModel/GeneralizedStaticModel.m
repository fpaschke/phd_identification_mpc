classdef GeneralizedStaticModel < idModels.StaticModel
    %GENERALIZEDSTATICMODEL Defines a simple static Model of kind 
    %   y = f(u,p) + e 
    %where e is assumed to be white noise with fixed covariancematrix Lambda_e (NoiseVariance).
    %A GeneralStaticModel can be of almost every form e. g. a polynomial or a exponential...
    %It is characterized by its Function handles "Fun" and its parametervectors "Parameters".
    %
    % obj = GeneralizedStaticModel(Nin,Nout,varargin)
    % obj [idModels.GeneralizedStaticModel]: GeneralizedStaticModel object
    % Nin [pos. int]: Number of inputs
    % Nout [pos. int]: Number of outputs
    % varargin: Name Value Pair Arguments: Any writable Property of IdModel/GeneralizedStaticModel to be set by Name value pairs                            
    %
    % LSQ-FIT for polynomial model with 2 inputs
    % hdl = @(u,p) deal(p(1)*u(:,1).^3 + p(2)*u(:,1).*u(:,2),[u(:,1).^3 u(:,1).*u(:,2)]);
    % u = 3*rand(1000,2);
    % [y,dy] = hdl(u,[3 2]);
    % y = y + 1*randn(1000,1) + 100;
    % GSM1 = idModels.GeneralizedStaticModel(2,1,'Fun',hdl,'Parameters',[3 1],'Name','GSM1','Free',[false true],'OutputOffset',50);
    % GSM1.identify(y,u,'CheckGradient',true,'EstimateOutputOffset',true);
    % GSM1.Parameters{1}; % should be [3 approx 2]
    % GSM1.showModel;
    %
    % ABSULUTE VALUE FIT FOR POLYNOMIAL WITH OUTLIERS
 	% hdl = @(u,p) deal(p(1)*u(:,1) + p(2)*u(:,1).^2,[u(:,1) u(:,1).^2]);
    % u = 3*rand(100,1);
    % [y,dy] = hdl(u,[2 3]);
    % y = y + 1*randn(length(u),1) + 100; 
    % y(50) = 300;
    % GSM2 = idModels.GeneralizedStaticModel(1,1,'Fun',hdl,'Parameters',[1 1],'Name','GSM2','Free',[true true],'OutputOffset',50);
    % GSM2.identify(y,u,'CheckGradients',true,'EstimateOutputOffset',true);
    % GSM2.showModel('Data',{y u});
    % GSM2.identify(y,u,'EstimateOutputOffset',true,'cost','abs','UseGradient',true);
    % GSM2.showModel('Data',{y u});
    
    
    properties
        Parameters % [ny x 1 cell of double vectors]: The Parameter Vectors of each output
        Free % [ny x 1 cell of double logicals]: The Parameter Vectors of each output
        Fun  %[ny x 1 cell of char]: The Functionhandles with two inputs describing the model f = @(u,p). u is considered to be a [N x nu] Matrix 
    end
    
    methods %public methods
        function obj = GeneralizedStaticModel(varargin) % Init from  Name Value Pairs
            if ~isempty(varargin) && isnumeric(varargin{1}) 
                Nin = varargin{1}; 
                varargin(1) = []; 
            else
                ix = find(strcmp(varargin,'InputName'));
                if ~isempty(ix) && ~ischar(varargin{ix+1})
                    Nin = length(varargin{ix+1});
                else
                    Nin = 1;
                end
            end
            if ~isempty(varargin) && isnumeric(varargin{1}) 
                Nout = varargin{1}; 
                varargin(1) = []; 
            else
                ix = find(strcmp(varargin,'OutputName'));
                if ~isempty(ix) && ~ischar(varargin{ix+1})
                    Nout = length(varargin{ix+1});
                else
                    Nout = 1;
                end
            end
            assert(isscalar(Nin) && Nin>0 && mod(Nin,1)==0,'"Nin" needs to be a pos int scalar!');
            assert(isscalar(Nout) && Nout>0 && mod(Nout,1)==0,'"Nout" needs to be a pos int scalar!');
            inname = arrayfun(@(i) ['u' num2str(i)],1:Nin,'UniformOutput',0)'; 
            outname = arrayfun(@(i) ['y' num2str(i)],1:Nout,'UniformOutput',0)';
            
            obj = obj@idModels.StaticModel('InputName',inname,'OutputName',outname);     
            obj = obj.set(varargin{:}); %Now set user supplied opt. name value pairs!
            if isempty(obj.Parameters)
                obj.Parameters = cell(obj.OutputDimension,1);
            end
           	if isempty(obj.Fun)
                obj.Fun = cell(obj.OutputDimension,1);
            end
        end
        
        %% Property Getters and Setters
        %Parameters
        function set.Parameters(obj,p)
            if obj.OutputDimension == 1 && isnumeric(p)
                assert(isnumeric(p) && all(~isnan(pi)),'The Parametervector needs to be a numeric value or a cell of numeric values!');
                obj.Parameters = {p};
            else
                assert(isvector(p) && length(p)==obj.OutputDimension,'The Parameter Vector needs to be a cell vector with as many cells as the model has outputs!');                
                assert(iscell(p) && all(cellfun(@(pi) all(~isnan(pi)) && isnumeric(pi) && (isvector(pi) || isempty(pi)),p)),'The Parameter Vector needs to be a cell with numeric vectors without NaN values!');
                obj.Parameters = p(:);
            end
            % Check weather allowed Parametervalue in case of subclass
            assert(obj.checkParameters,'The Value youre trying to assign to parameters does not match the requirements! Check "checkParameters" method of your model class!');
            
            % Check length of Free
            for no = 1:obj.OutputDimension
                if isempty(obj.Free) || length(obj.Parameters{no})~=length(obj.Free{no}) 
                    free{no} = true(size(obj.Parameters{no}));
                else
                    free{no} = obj.Free{no};
                end
            end
            obj.Free = free;
        end
        function p = get.Parameters(obj)
            p = obj.Parameters;
        end
       	%Free
        function set.Free(obj,p)
            if obj.OutputDimension == 1 && ~iscell(p)
                p = {p};
            end
            p = cellfun(@(pi) logical(pi),p,'UniformOutput',false);
            assert(length(p) == length(obj.Parameters) && all(arrayfun(@(no) islogical(p{no}) && length(p{no}) == length(obj.Parameters{no}),...
                1:obj.OutputDimension)),'The Free vector needs to be a cell of logical vectors with same Dimensions as obj.Parameters!');  
            obj.Free = p(:);
        end
        function p = get.Free(obj)
            p = obj.Free;
        end
      	%Fun
        function set.Fun(obj,fs)
            if isa(fs,'function_handle') fs = func2str(fs); end
            if ischar(fs) fs = {fs}; end
            if all(cellfun(@isempty,fs))
                obj.Fun = fs(:);
                return;
            end
            assert(isvector(fs) && length(fs)==obj.OutputDimension,'The Value to be set needs to be a cell vector of string function handles with as many cells as the model has outputs!');                
            f = cellfun(@str2func,fs,'UniformOutput',0);
            assert(iscell(f) && all(cellfun(@(fi) isa(fi,'function_handle') && nargin(fi)>=2,f)),'All function handles need at least two inputs!');
            obj.Fun = fs(:);
        end
        function f = get.Fun(obj)
            f = obj.Fun;
        end
        
        %% OTHER
        % Abstract implementations
        opt = identifyModelOptions(obj,varargin)
        e = identifyModel(obj,y,u,opt)
        y = simulateModel(obj,u,opt)
    end %public methods
    
    methods (Access = private)
      	function setPvec(obj,P)
            ixs = 1;
            for no = 1:obj.OutputDimension
                ixe = ixs + length(obj.Parameters{no}(obj.Free{no}))-1;
                obj.Parameters{no}(obj.Free{no}) = P(ixs:ixe);
                ixs = ixe + 1;
            end
        end 
        
        function [freeP, free, P] = getPvec(obj)
            P = []; free = logical([]);
            for no = 1:obj.OutputDimension
                P = [P obj.Parameters{no}];
                free = [free obj.Free{no}];
            end
            P = P';
            freeP = P(free);
        end 
    end
    
 	%% PROTECTED METHS
    methods (Access = protected)

    end
  	%% PROTECTED ABSTRACT METH
    methods (Access = protected)
        function ok = checkParameters(obj)
            ok = true; %No checks can be made here
        end
    end
    
  	%% Static Methdos
    methods (Static = true)
        function obj = loadobj(str)
            nvp = util.struct2namevaluepairs(str);
            % Make Sure to set Parameters first! (Before Free!)
            ix = find(strcmp(nvp,'Parameters'));
            val = nvp{ix+1};
            nvp([ix ix+1]) = [];
            nvp = {'Parameters' val nvp{:}};
            obj = idModels.GeneralizedStaticModel(length(str.InputName),length(str.OutputName),nvp{:});
        end
    end
end
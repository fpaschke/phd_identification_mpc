classdef (Abstract) IdModel < handle
    %IDMODEL Defines a abstract Model from which actual model classes will be derived. 
    
    properties
        Name = 'unnamed'    % [string]: Name of a Model
        Info = [];          % [struct]: Infostruct contains Information about the model
        InputName = {};     % [nu x 1 cell of strings]: InputNames of Model
        InputUnit = {};     % [ny x 1 cell of strings]: InputUnits of Model 
        InputMin = [];      % [nu x 1 double]: Minimum of Inputsignals 
        InputMax = [];      % [nu x 1 double]: Maximum of Inputsignals
        OutputName = {};    % [ny x 1 cell of strings]: OutputNames of Model 
        OutputUnit = {};    % [ny x 1 cell of strings]: OutputUnits of Model 
        OutputMin = [];     % [ny x 1 double]: Minimum of Outputsignals
        OutputMax = [];     % [ny x 1 double]: Maximum of Outputsignals
        OutputOffset = [];  % [ny x 1 double]: of outputoffsets y0 y = y_ + y0 [Def.: 0]
       	NoiseVariance = []; % [ny x ny double]: Covariance matrx of outputnoise (Innovations)
    end
    
    properties (Dependent)
        InputDimension      % [pos int scalar]: Number of Model Inputs 
        OutputDimension     % [pos int scalar]: Number of Model Outputs
    end
    
    %% PUBLIC METHODS
    methods
        %% Constructor
        function obj = IdModel(varargin)
            % No desccriprion of constructor available since IdModel cannot be instantiated.
            % Make Sure to set InputName and OutputName first!
            ix = find(strcmpi(varargin,'InputName'));
            if ~isempty(ix)
                obj = set(obj,'InputName',varargin{ix+1});
            end
            ix = find(strcmpi(varargin,'OutputName'));
        	obj = set(obj,'OutputName',varargin{ix+1});

            obj = set(obj,varargin{:});
          	if isempty(obj.NoiseVariance); obj.NoiseVariance = NaN(obj.OutputDimension); end
            if isempty(obj.InputMin); obj.InputMin = -Inf(obj.InputDimension,1); end
            if isempty(obj.InputMax); obj.InputMax = Inf(obj.InputDimension,1); end
            if isempty(obj.OutputMin); obj.OutputMin = -Inf(obj.OutputDimension,1); end
            if isempty(obj.OutputMax); obj.OutputMax = Inf(obj.OutputDimension,1); end
            if isempty(obj.OutputOffset); obj.OutputOffset = zeros(obj.OutputDimension,1); end
            if isempty(obj.InputMin); obj.InputMin = -Inf(obj.InputDimension,1); end
            if isempty(obj.InputUnit); obj.InputUnit = arrayfun(@(x) '',1:obj.InputDimension,'UniformOutput',false); end
            if isempty(obj.OutputUnit); obj.OutputUnit = arrayfun(@(x) '',1:obj.OutputDimension,'UniformOutput',false); end
        end
        
        %% Property Getters and Setters
        % Name
        function set.Name(obj,name)
            assert(ischar(name),'The Name needs to be a string!');
            obj.Name = name;
        end
        function name = get.Name(obj)
            name = obj.Name;
        end
       	% InputName
        function set.InputName(obj,lab)
            assert(ischar(lab) || iscellstr(lab),'InputName needs to be a string or a cell vector of strings!');
            if ischar(lab); lab = {lab}; end
            assert(length(lab) == length(unique(lab)),'All Input Names need to be different!');
            assert(isvector(lab) || isempty(lab),'InputName needs to be a vector!');
            assert(length(lab)==length(obj.InputName) || isempty(obj.InputName),'Number of Inputs is not allowed to change!');
            obj.InputName = lab(:);
        end
        function lab = get.InputName(obj)
            lab = obj.InputName;
        end
        % OutputName
        function set.OutputName(obj,lab)
            assert(ischar(lab) || iscellstr(lab),'OutputName needs to be a string or a cell of strings!');
            if ischar(lab); lab = {lab}; end
            assert(length(lab) == length(unique(lab)),'All Output Names need to be different!');
            assert(isvector(lab),'OutputName needs to be a vector!');
            assert(length(lab)==length(obj.OutputName) || isempty(obj.OutputName),'Number of Inputs is not allowed to change!');
            obj.OutputName = lab(:);
        end
        function lab = get.OutputName(obj)
            lab = obj.OutputName;
        end
    	% InputUnit
        function set.InputUnit(obj,un)
            assert(ischar(un) || iscellstr(un),'InputUnit needs to be a string or a cell vector of strings!');
            if ischar(un); un = {un}; end
            assert(isvector(un) || isempty(un),'InputUnit needs to be a vector!');
            assert(length(un)==obj.InputDimension,'Number of input units needs to be the same as the number of Inputs!');
            obj.InputUnit = un(:);
        end
        function un = get.InputUnit(obj)
            un = obj.InputUnit;
        end
        % OutputUnit
        function set.OutputUnit(obj,un)
            assert(ischar(un) || iscellstr(un),'OutputUnit needs to be a string or a cell vector of strings!');
            if ischar(un); un = {un}; end
            assert(isvector(un) || isempty(un),'OutputUnit needs to be a vector!');
            assert(length(un)==obj.OutputDimension,'Number of input names needs to be the same as the number of Outputs!');
            obj.OutputUnit = un(:);
        end
        function un = get.OutputUnit(obj)
            un = obj.OutputUnit;
        end
        % InputMin
        function set.InputMin(obj,val)
            assert(isvector(val) || isempty(val),'The value to be set needs to be a vector!');
            val = val(:);
            assert(length(val) == obj.InputDimension,'The value to be set needs to be a vector of obj.InputDimension!');
            assert(isnumeric(val) && all(isreal(val)) && all(~isnan(val)),'All values need to be numeric real values!');
            if ~isempty(obj.InputMax); assert(all(val<=obj.InputMax),'InputMin should be <= InputMax!'); end
            obj.InputMin = val;
        end
        function val = get.InputMin(obj)
            val = obj.InputMin;
        end
        % InputMax
        function set.InputMax(obj,val)
            assert(isvector(val) || isempty(val),'The value to be set needs to be a vector!');
            val = val(:);
            assert(length(val) == obj.InputDimension,'The value to be set needs to be a vector of obj.InputDimension!');
            assert(isnumeric(val) && all(isreal(val)) && all(~isnan(val)),'All values need to be numeric real values!');
            if ~isempty(obj.InputMin); assert(all(obj.InputMin<=val),'InputMin should be <= InputMax!'); end
            obj.InputMax = val;
        end
        function val = get.InputMax(obj)
            val = obj.InputMax;
        end
        % OutputMin
        function set.OutputMin(obj,val)
            assert(isvector(val),'The value to be set needs to be a vector!');
            val = val(:);
            assert(length(val) == obj.OutputDimension,'The value to be set needs to be a vector of obj.OutputDimension!');
            assert(isnumeric(val) && all(isreal(val)) && all(~isnan(val)),'All values need to be numeric real values!');
            if ~isempty(obj.OutputMax); assert(all(val<=obj.OutputMax),'OutputMin should be <= OutputMax!'); end
            obj.OutputMin = val;
        end
        function val = get.OutputMin(obj)
            val = obj.OutputMin;
        end
        % OutputMax
        function set.OutputMax(obj,val)
            assert(isvector(val),'The value to be set needs to be a vector!');
            val = val(:);
            assert(length(val) == obj.OutputDimension,'The value to be set needs to be a vector of obj.OutputDimension!');
            assert(isnumeric(val) && all(isreal(val)) && all(~isnan(val)),'All values need to be numeric real values!');
            if ~isempty(obj.OutputMin); assert(all(obj.OutputMin<=val),'OutputMin should be <= OutputMax!'); end
            obj.OutputMax = val;
        end
        function val = get.OutputMax(obj)
            val = obj.OutputMax;
        end
        %NoiseVariance
        function set.NoiseVariance(obj,sigma2)
            if isvector(sigma2)
                sigma2 = diag(sigma2);
            end
            assert(all(all(isreal(sigma2) & isequaln(sigma2,sigma2'))) & (size(sigma2,1) == obj.OutputDimension), ...
                'NoiseVariance needs to be symetric matrix with as many rows/columns as the model has outputs!');
            obj.NoiseVariance = sigma2;
            if isa(obj,'idModels.NsfPolyModel') %Invoke rebuild of Ss object. 
                obj.updateSs(); 
            end
        end
        function sigma2 = get.NoiseVariance(obj)
            sigma2 = obj.NoiseVariance;
        end
        % Info
        function set.Info(obj,I)
            assert(isstruct(I) || isempty(I),'Info needs to be a structure!');
            obj.Info = I;
        end
        function I = get.Info(obj)
            I = obj.Info;
        end
        %OutputOffset
        function set.OutputOffset(obj,val)
            if ~isempty(val)
                assert(isvector(val) && length(val) == obj.OutputDimension,'OutputOffset needs to have as many values as the model has outputs!');
                assert(all(isreal(val)) && all(isnumeric(val)) && all(~isinf(val)) && all(~isnan(val)),'OutputOffset needs to be real numric values or empty!');
                obj.OutputOffset = val(:);
            else
                obj.OutputOffset = zeros(obj.OutputDimension,1);
            end
        end
        function y0 = get.OutputOffset(obj)
            y0 = obj.OutputOffset;
        end
        %% Dependent Props
      	% InputDimension
        function m = get.InputDimension(obj)
            m = length(obj.InputName);
        end
        % OutputDimension
        function p = get.OutputDimension(obj)
            p = length(obj.OutputName);
        end
        
        %% UTIL
    	function [y,u] = getRawData(obj,y,u)
            % if y is a tsc the extract raw I/O data 
            if isa(y,'signals.TSignalCollection')
                [y,u,~] = obj.getRawDataFromTsc(y,true);
            end
            
            % bring y and u to cells
            if ~iscell(y); y = {y}; end
            if size(u,1)==0; u = cellfun(@(yi) zeros(length(yi),0),y,'UniformOutput',false); end
            %if isempty(u); u = cellfun(@(yi) zeros(length(yi),0),y,'UniformOutput',false); end
            if ~iscell(u); u = {u}; end
            y = y(:); u = u(:);
            assert(all(obj.InputDimension==cellfun(@(ui) size(ui,2),u)),'The number of supplied input signals does not match the number of model inputs!');
            if ~isempty(y{1}) && length(y) == 1 % if y = [] (y isnt necessary, eg for simulation) 
                assert(all(obj.OutputDimension==cellfun(@(yi) size(yi,2),y)),'The number of supplied output signals does not match the number of model outputs!');
            end
        end
    end %public methods
    
    %% Hidden and protected methods
    methods (Hidden = true, Access = protected)   
        function obj = set(obj,varargin)
            %SET Sets Properties of obj with name value pairs
            %
            % obj = set(obj,varargin)
            % obj [idModel]: a idModel object
            % varargin: Valid Name Value Pairs of obj to be set, or structure with valid name value pairs.
            if length(varargin) == 1 % Case struct
                if isstruct(varargin{1})
                    fn = fieldnames(varargin{1});
                    for np = 1:length(fn)
                        obj.(fn{np}) = varargin{1}.(fn{np});
                    end
                else
                    error('Set works with name value pairs or wioth a structure with properties to be set!');
                end
            else    
                assert(mod(length(varargin),2)==0,'The Number of input arguments needs to be even (Name Property values only)!');
                assert(iscellstr(varargin(1:2:end)),'All Inputs need to be Name Property Values!');
                assert(length(unique(varargin(1:2:end)))==length(varargin(1:2:end)),'Each Name Value pair can only be assigned once!');

                % Set all available Properties
                for np = 1:2:length(varargin)
                    obj.(varargin{np}) = varargin{np+1};
                end
            end
        end
        
        function checkData(obj,y,u)
            %CHECKDATA Checks i/o data for Model Idenfication or Simulation.
            %
            % checkInputs(obj,y,u,varargin)
            % obj [IdModel]: Any Valid BModel Object
            % y [(cell of) N x ny Matrix]: Matrix of Measured Output Data
            % u [opt. (cell of) N x nu Matrix]: Matrix of Measured Input Data

            if ~iscell(y) && ~isempty(y); y = {y}; end
            if nargin>2 && ~iscell(u) && ~isempty(u); u = {u}; end

            if nargin>2 && ~isempty(y) && ~isempty(u) && any(size(y) ~= size(u))
                error('Sizes of u an y need to be equal!');
            end

            if ~isempty(y)
                for ns = 1:length(y)
                        assert(ismatrix(y{ns}),'y need to be matrix!')
                        assert(isequal(size(y{ns},2),obj.OutputDimension),'Dimensions of Outputdata need to match OutputDimension of Model!');
                        assert(isnumeric(y{ns}) && all(all(isreal(y{ns}))),'y are allowed to contain only real numeric values (non NaN)!');
                        assert(all(all(~isnan(y{ns}))),['OutputDataset ' num2str(ns) ' contains NaNs']);
                end
            end
            if nargin>2 && ~isempty(u)
                for ns = 1:length(u)
                    assert(ismatrix(u{ns}),'u need to be matrix!')
                    assert(isequal(size(u{ns},2),obj.InputDimension),'Dimensions of Inputdata needs to match InputDimension of Model!');
                    assert(isnumeric(u{ns}) && all(all(isreal(u{ns}))),'u are allowed to contain only real numeric values (non NaN)!');
                    assert(all(all(~isnan(u{ns}))),['InputDataset ' num2str(ns) ' contains NaNs']);
                end
            end
        end
        
        function [y, u, tsc] = getRawDataFromTsc(obj,tsc,eco)
            %Function extracts selected data form TSignalCollection for identification or simulation.
            %NOTE: The Data will be selected by Identifiers! This means that the function that obj.OutputName
            %and obj.InputName are assumed to be the identifiers of the signals contained in tsc.
            %
            % [y, u, tsc] = getRawDataFromTsc(obj,tsc,eco)
            % obj [IdModel]: IDModel object for which simulation/identification data shall be extracted
            % tsc [TSignalCollection]: Collection containing the Signals necessary for simulation/identification.
            % eco [opt. 0/1]:   Eco Mode: This flag only influences the output tsc. 
            %               	If 1 then only the IO Signals of obj will be resampled/synchronized and saved to tsc. (Default: 1).
            %                   All other signals in tsc will be kept unchanged. 0 will resample synchronize and select all signals in tsc.                	
            % [y (cell of) Ns x ny double]: Extracted outputdata. 
            % [u (cell of) Ns x nu double]: Extracted inputdata. 
            % tsc [TSignalCollection]: Resampled/Synchronized and Selected data with IO Data of tsc only.

            % subset data
            nu = obj.InputDimension;
            ny = obj.OutputDimension;
            io_sigs = [obj.InputName' obj.OutputName'];
            if nargin > 2 && ~eco
                tsc_io = tsc;
            else
                eco = 1;
                tsc_io = tsc.getSignalsFromIdentifiers(io_sigs);
            end

            % check if synchronized
            if isa(obj,'idModels.StaticModel') && tsc_io.Length > 1  
                warning([   'The Input and Outputsignals of the Model "%s" did not have the same timebase.' ...
                            'The Signals will thus be synchronized!'],obj.Name);
                tsc_io = tsc_io.synchronize();
            elseif isa(obj,'idModels.DynamicModel') % For dynamic models a resampling to Sampling Time of obj needs to be done!
                % sampling times and timeunits of the signals
                Tu  = cellfun(@(n) tsc_io.(n).TimeInfo.Unit,tsc_io.getSignalNames,'UniformOutput',false);
                if length(unique(Tu))~=1
                    tsc_io = tsc_io.convertTimeUnits(Tu{1});
                end
                Sd  = cellfun(@(n) tsc_io.(n).TimeInfo.StartDate,tsc_io.getSignalNames,'UniformOutput',false);
                if ~all(cellfun(@isempty,Sd))
                    if length(unique(datenum(Sd)))~=1
                        tsc_io = tsc_io.convertStartDates(Sd{1});
                    end
                end

                % convert samplingtime of the model into timeunit of tsc_io
                Ts = cellfun(@(n) tsc_io.(n).SamplingTime,tsc_io.getSignalNames);
                Tsmod = obj.Ts*time.getTimeUnitConvFactor(obj.TimeUnit)/time.getTimeUnitConvFactor(Tu{1}); 
                if any(abs(Ts-Tsmod)>10*eps) || any(isnan(Ts)) || tsc_io.Length > 1  
                    warning([   'At least one of the Input-/Outputsignals in the TSignalCollection did not have the samplingtime of the Model "%s". ' ...
                                'All Signals will thus be resampled to %d %s!'],obj.Name,obj.Ts,obj.TimeUnit);

                    Ti = cell2mat(cellfun(@(n) [tsc_io.(n).Time(1) tsc_io.(n).Time(end)]',tsc_io.getSignalNames','UniformOutput',false));
                    tsc_io = tsc_io.resample(max(Ti(1,:)):Tsmod:min(Ti(2,:))); % 28.4.2020: problematisch z B wenn signale nur aus einem Sample bestehen
                    % tsc_io = tsc_io.resample(min(Ti(1,:)):Tsmod:max(Ti(2,:)));
                end
            end

            % Due to resampling/synchronization the selection vectors are not necessarilly the same
            % we thus compute an new/global selection by ORing all the selection vectors
            names = tsc_io.getSignalNames;
            for ns = 1:length(names)
               S(:,ns) = tsc_io.(names{ns}).Selection; 
            end
            if nnz(diff(S,1,2))~=0
                warning([   'Unequal Selection Vectors appeared after resampling/synchronization! ' ... 
                            'A new Selection vector for all resampled/synchronized Signals will ' ...
                            'be computed and set to all TSignals using "OR" Operation!']);
            end
            SOR = any(S,2); %columnwise or
            for ns = 1:length(names)
               tsc_io.(names{ns}).Selection = SOR; 
            end

            % now get the Data
            yu = tsc_io.getDataFromSelection(1,'Identifiers',[obj.OutputName; obj.InputName]);
            if nu>0
                [y, u] = arrayfun(@(ns) deal(yu{ns}(:,1:ny),yu{ns}(:,ny+1:end)),(1:length(yu))','UniformOutput',false);
            else
                y = yu;
                u = cell(length(y),1);
            end

            if eco % write tsc_io to tsc
                no_io_signals = tsc.getIdentifierNames;
                for k = 1:length(io_sigs)
                    no_io_signals = no_io_signals(~strcmp(no_io_signals,io_sigs{k}));
                end
                tsc = signals.TSignalCollection({tsc_io tsc(no_io_signals)});
            else % output
                tsc = tsc_io;
            end
        end
        
        function lb = getInLab(obj,num)
            lb = obj.InputName{num};
            if ~isempty(obj.InputUnit{num})
                lb = [lb ' [' obj.InputUnit{num} ']'];
            end
            lb = util.convToPlotLabel(lb,'latex');
        end

        function lb = getOutLab(obj,num,nm)
            if nargin < 3
                lb = obj.OutputName{num};
            else
                lb = nm;
            end     
            if ~isempty(obj.OutputUnit{num})
                lb = [lb ' [' obj.OutputUnit{num} ']'];
            end
            lb = util.convToPlotLabel(lb,'latex');
        end
    end
    
    
  	%% Abstract protected Methods: To be implemented in subclasses
    methods (Abstract = true, Hidden = true)
       y = simulateModel(obj,y,u,varargin);
       e = identifyModel(obj,y,u,varargin);
       opt = identifyModelOptions(obj);
    end
    
    %% Static Methods
    methods (Static = true)
        Lambda = calcCovariance(e,Npar,SubMean)  
    end
    
    %% Abstract static Methods
  	methods (Static = true, Abstract = true)
        loadobj(str)
    end
end
function e = identify(obj,y,varargin)
%IDENTIFY Identify a model with given IO Data y and u.
%
% Internal functionality:
% 1.)   Creates an appropriate option struct for identification of obj. See identifyModelOptions.
% 2.)   Prepare and check raw data for identification. 
% 3.)   Remove OutputOffset from y data.
% 4.)   Call the identifyModel of specific subclass to invoke parameter identification 
% 5.)   Calculate, print and write some statistics to obj.Info struct
% The specific identifyModel shall do the following:
% 1.)   Estimate parameters of model.
% 2.)   Write the parameters to the model object
%
% e = identify(obj,y,u,varargin)
%   obj [BModel]: Any subclass of BModel with an identifyModel method.
%   y [(cell of) N x ny double Matrix]: Raw Output Data for Model identification
%   u [(cell of) N x nu double Matrix]: Raw Input Data for Model identification
%   e [(cell of) N x ny double Matrix]: Residuals resp. 1-step ahead prediction errors
%   varargin: One of the following
%       opt [opt. Structure]: Structure Specifying the options to use for estimation -> see identifyOptions.m of the model to be identified. 
%   OR:
%       Valid Name Value Pairs: see identifyOptions.m for valid Options of the model to be identified. 

%assert(obj.OutputDimension == 1,'Identification for MultiOutput Systems is unsupported yet!');

%% INIT
% GET u
if ~isempty(varargin) && (isnumeric(varargin{1}) || iscell(varargin{1}))
    u = varargin{1}; varargin(1) =[];
else
    u = [];
end

% GET OPTIONS 
if length(varargin) >= 2 || isempty(varargin)
    opt = obj.identifyOptions(varargin{:});
elseif isstruct(varargin{1}) && length(varargin) == 1
    opt = util.struct2namevaluepairs(varargin{1});
    opt = obj.identifyOptions(opt{:});
else
   error('Invalid Input Argument!')
end
            
% GET RAW DATA       
if iscell(y); cellin = true; else; cellin = false; end
[y , u] = getRawData(obj,y,u);
obj.checkData(y,u);

% Set Minima and Maxima
obj.OutputMin = min(cell2mat(y));
obj.OutputMax = max(cell2mat(y));
obj.InputMin = min(cell2mat(u));
obj.InputMax = max(cell2mat(u));

%% DO IDENTIFICATION AND CALC COVARIANCE MATRIX
e = obj.identifyModel(y,u,opt);
%obj.NoiseVariance = obj.calcCovariance(e,obj.Info.NumEstParameters);

%% SET INFO AND PRINT STATS
e_mat = cell2mat(e);
ix_dat = all(~isnan(e_mat),2);
Info = obj.Info;
Info.EstimatedOn = datestr(now);
Info.OptionsUsed = opt;
Info.SamplesUsed = sum(length(cell2mat(e)));
Info.SamplesUsedCost = sum(ix_dat);
Info.MeanE = mean(e_mat(ix_dat,:),1)';
Info.Mse = mean(e_mat(ix_dat,:).^2,1)';
Info.Mae = mean(abs(e_mat(ix_dat,:)),1)';
obj.Info = Info;
fprintf('Performed parameter estimation for Model "%s" using %i Samples: \n',obj.Name,obj.Info.SamplesUsed);
fprintf('Mean of Noise: \n');
disp(obj.Info.MeanE);
fprintf('NoiseVariance: \n');
disp(obj.NoiseVariance);

if ~cellin; e = e{1}; end
end


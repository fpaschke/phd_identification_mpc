function opt = identifyOptions(obj,varargin)
%IDENTIFYOPTIONS Create Optionstructure for model identification. 
%
% opt = getIdentifyOptions(obj,varargin)
% obj [idModels.IdModel]: A instance of any subclass of IdModel
% varargin: Name Value Pairs defining options for identification. 
%   'EstimateOutputOffset' [logical]: If true the OutputOffset will be treated as an free parameter. Default: false.

p = inputParser(); 
p.KeepUnmatched = true;
% Define Options that are valid for identification of all subclasses of IdModel here
addParameter(p,'EstimateOutputOffset',false,@(x) all(x==0 | x==1) && length(x) == obj.OutputDimension);
parse(p,varargin{:});
opt = p.Results;

opt2add = obj.identifyModelOptions(p.Unmatched);
% Append optional args 
f = fieldnames(opt2add);
for i = 1:length(f)
    opt.(f{i}) = opt2add.(f{i});
end
end
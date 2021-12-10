function [e,yp,stats] = calcResiduals(obj,y,u,varargin)
%CALCRESIDUALS Calculates residuls of the Model by comparing simulated model response to measured data.
%   
% [e,yp,misc] = calcResiduals(obj,y,u)
% obj [IdModel]: Any valid model for which simulate routine is available.
% y [(cell of) N x ny double]: Measured OutputData
% u [(cell of) N x nu double]: Measured InputData
% varargin: optional name value pairs to passed to sumulate (static Model) or calcPredictions (NsfPolyModel)
% e [(cell of) N x ny x Hp double]: Residuals aka prediction errors. 
% yp [(cell of) nu x N double]: Measured InputData
% stats [struct]: Structure containing statistics such as MSE, MAE, STD, etc. of residuals e

%% INIT
obj.checkData(y,u);
if ~iscell(u); cellin = false; else; cellin = true; end

%% SIMULATE SYSTEM
if isa(obj,'idModels.StaticModel')
    yp = obj.simulate(u,varargin{:}); 
elseif isa(obj,'idModels.NsfPolyModel') % call vectorized implementation of simulate.
    yp = obj.calcPredictions(y,u,varargin{:});
else
    error('Not Implemented yet')
end

%% CALC RESIDUALS
if ~cellin; yp = {yp}; end
if ~iscell(y); y = {y}; end

Hp = size(yp,3);
Nsets = length(y); 
e = cell(Nsets,1);
for ns = 1:Nsets
    ym = repmat(y{ns},[1 1 Hp]);
    e{ns} = yp{ns} - ym;
end

%% CALC STATS
if nargout > 2
    E = cell2mat(e);
    allErr = E(~isnan(E(:,end)),:,:);
    stats.Npred = size(allErr,1);
    stats.Mean = squeeze(mean(allErr))';
    stats.Mae = squeeze(mean(abs(allErr)))';
    stats.Mse = squeeze(mean(allErr.^2))';
    stats.Std = squeeze(std(allErr))';
    E = E(:); E = E(~isnan(E));
    stats.Rmmse = sqrt(mean(E.^2)); % RMMSE ERROR

    if obj.OutputDimension == 1
        stats.Mean = stats.Mean'; stats.Std = stats.Std'; stats.Mae = stats.Mae'; stats.Mse = stats.Mse';
    end
end
%% OUTPUT
if ~cellin; e = e{1}; yp = yp{1}; end
end
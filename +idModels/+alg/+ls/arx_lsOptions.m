function opt = arx_lsOptions(varargin)
%Generates options struct for arx_ls
%
% opt = arx_lsOptions(varargin)
% varargin: Name Value Pairs
%   'Prefilter' [ny cell of 1x2 cell]: cell array with two vectors defining a Prefilter/WeightingFilter. first cell determines numerator 
%                           coefficients and 2nd contains denum coefficients of the filter. Both need to be monic polynomials
%                           Def. [] (no prefiltering)
%   'Instruments' [ny cell or double matrix]:   Instrumental Variables to be used. Shall have the same format / dimension as
%                                               output data. Def. [] (no instruments will be used)

p = inputParser();
addOptional(p,'Instruments',[]);
valFun = @(x) isempty(x) || isreal(x) || (iscell(x) && length(x(:)) == 2 && x{1}(1)==1 && x{2}(1)==1 && all(isnumeric([x{1} x{2}])));
addParameter(p,'Prefilter',[],@(x) isempty(x) || all(cellfun(@(pi) valFun(pi),x)));
parse(p,varargin{:});
opt = p.Results;
end
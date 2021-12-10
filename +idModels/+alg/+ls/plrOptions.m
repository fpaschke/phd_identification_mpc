function opt = plrOptions(varargin)
%Generates options struct for plr.
%
% opt = plrOptions(varargin)
% varargin: Name Value Pairs
%   'Prefilter' [ny cell of 1x2 cell]: cell array with two vectors defining a Prefilter/WeightingFilter. first cell determines numerator 
%                           coefficients and 2nd contains denum coefficients of the filter. Both need to be monic polynomials
%                           Def. [] (no prefilter)
%   'Plot' [log]: Creates plot of polynomial convergence if true. Def. false
%   'Iter' [log]: Creates plot of polynomial convergence if true. Def. 50
%   'Init' ['ls' or 'iv']: Determines the initialization step. 'ls' will fit an arx model via Least squares first, while 'iv' will use instrumental variable method. Def. 'iv'

p = inputParser();
addParameter(p,'Plot',false,@(x) x==1 || x==0); 
addParameter(p,'Iter',20,@(x) x>=1 && mod(x,1) == 0); 
addParameter(p,'Init','ls',@(x) any(strcmpi(x,{'iv' 'ls'}))); 
valFun = @(x) isempty(x) || isreal(x) || (iscell(x) && length(x(:)) == 2 && x{1}(1)==1 && x{2}(1)==1 && all(isnumeric([x{1} x{2}])));
addParameter(p,'Prefilter',[],@(x) isempty(x) || all(cellfun(@(pi) valFun(pi),x)));
parse(p,varargin{:});
opt = p.Results;
end
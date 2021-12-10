function opt = identifyModelOptions(obj,varargin)
%IDENTIFYMODELOPTIONS Generates options struct for identify function of GeneralizedStaticModel.

p = inputParser; p.KeepUnmatched = true;
% Add Parameters here
addParameter(p,'Display','iter');
addParameter(p,'Solver','matlab_lsqnonlin',@ischar);
addParameter(p,'UseGradient',true,@(x) x==0 || x==1);
addParameter(p,'CheckGradient',false,@(x) x==0 || x==1);
addParameter(p,'Cost','quad',@(x) any(strcmpi(x,{'quad' 'abs'})));

parse(p,varargin{:});
opt = p.Results;
opt.Unmatched = p.Unmatched;
end
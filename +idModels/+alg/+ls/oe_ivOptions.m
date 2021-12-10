function opt = oe_ivOptions(varargin)
%OE_IVOPTIONS Creates Options structure for oe_iv function.

p = inputParser();
valFun = @(x) isempty(x) || isreal(x) || (iscell(x) && length(x(:)) == 2 && x{1}(1)==1 && x{2}(1)==1 && all(isnumeric([x{1} x{2}])));
addParameter(p,'Iter',20,@(x) x>0 && mod(x,1) == 0); 
addParameter(p,'Plot',false,@(x) x==1 || x==0); 
addParameter(p,'EstIc',true,@(x) x==1 || x==0); 
addParameter(p,'Ls_init_samples','auto',@(x) (isnumeric(x) && (isinf(x) || mod(x,1)==0)) || strcmpi(x,'auto'));
addParameter(p,'StabilizationFactor',1.1,@(x) x>1);
addParameter(p,'Prefilter',[],@(x) isempty(x) || all(cellfun(@(pi) valFun(pi),x)));
parse(p,varargin{:});
opt = p.Results;
end
function [F,dFdP] = fun_logdThcv(u,p,y)
% Logarithmic function multiplied by dT
%
%       ln(exp(p1)*(u1) + 1)
% F = ---------------------------------------  ( [fun_softplus(u(:,end)) + p(end)] - Troom)^p(2) 
%               ln(exp(p1) + 1) 
%

% p(1): Scaling for valves
% p(2): Heizkörperexp
% p(3): Anstieg von hcv -> fix
% p(4): Scaling von hcv
% p(5): Shift von hcv
% p(6): offset hcv

N = size(u,1);
if nargin < 3
    h = u(:,1:end-2);
    [Tsup ,dTsup] = idModels.func.fun_softplus(u(:,end-1),p(3:end-1));
    Tr = u(:,end-1);
else
    h = u(:,1:end-1);
    [Tsup ,dTsup] = idModels.func.fun_softplus(u(:,end),p(3:end-1));
    Tr = y;
end
    
[f1, df1] = idModels.func.fun_log(h,p(1));
[f2, df2] = idModels.func.fun_dT(Tsup + p(end),p(2),Tr);
F = f1.*f2;                        
dFdP = [f2.*df1 f1.*df2 f1.*[dTsup ones(N,1)].*p(2).*idModels.func.fun_dT(Tsup + p(end),p(2)-1,Tr)];
if any(~isreal(F)) || any(isnan(dFdP(:)))
    stop = 1;
end

end
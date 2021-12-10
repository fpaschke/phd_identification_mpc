function [F,dFdP] = fun_polydT(u,p,y)
% Logarithmic function multiplied by dT
%
%       ln(exp(p1)*(u1) + 1)
% F = ---------------------------------------  (Tsup - Troom)^p(2) 
%               ln(exp(p1) + 1) 
%

if nargin < 3
    h = u(:,1:end-2);
    Tsup = u(:,end-1);
    Tr = u(:,end-1);
else
    h = u(:,1:end-1);
    Tsup = u(:,end);
    Tr = y;
end
    
[f1, df1] = idModels.func.fun_poly(h,p(1:end-1));
[f2, df2] = idModels.func.fun_dT(Tsup,p(end),Tr);
F = f1.*f2;
dFdP = [f2.*df1 f1.*df2];
if any(~isreal(F)) || any(isnan(dFdP(:)))
    stop = 1;
end

end
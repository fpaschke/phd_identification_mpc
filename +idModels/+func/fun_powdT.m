function [F,dFdP] = fun_powdT(u,p,y)
% Power function multiplied by dT
%
% narargin == 3:
%
% F = mean(u(:,1:end-1),2).^p(1) (u(:,end) - y)^p(2) 
%
% narargin == 2:
%
% F = mean(u(:,1:end-2),2).^p(1) (u(:,end) - u(:,end-1))^p(2) 
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
    
[f1, df1] = idModels.func.fun_pow(h,p(1));
[f2, df2] = idModels.func.fun_dT(Tsup,p(2:end),Tr);
F = f1.*f2;
dFdP = [f2.*df1 f1.*df2];
% if any(isnan(F)) || any(~isreal(F)) || any(isnan(dFdP(:))) || any(~isreal(dFdP(:)))
%    stop = 1; 
% end
end
function [F,dFdP] = fun_logdT(u,p,y)
% Logarithmic function multiplied by dT
%
% narargin == 3:
%
% F = ln(exp(p1)*mean(u(:,1:end-1),2) + 1) * (u(:,end) - y)^p(2) 
%
% narargin == 2:
%
% F = ln(exp(p1)*mean(u(:,1:end-2),2) + 1) *  (u(:,end) - u(:,end-1))^p(2) 
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
    
[f1, df1] = idModels.func.fun_log(h,p(1));
[f2, df2] = idModels.func.fun_dT(Tsup,p(2:end),Tr);
F = f1.*f2;
dFdP = [f2.*df1 f1.*df2];
% if any(isnan(F)) || any(~isreal(F)) || any(isnan(dFdP(:))) || any(~isreal(dFdP(:)))
%    stop = 1; 
% end
end

% u2 = -10:0.1:10; u1 = linspace(0,1,length(u2)); dp = 1e-8; p = [2 1]; T = zeros(length(u1),1)
% [y,dy] = idModels.func.fun_logdT([u1' u2'],p,T);
% dy_1 = (idModels.func.fun_logdT([u1' u2'],p+dp*[1 0],T) - idModels.func.fun_logdT([u1' u2'],p-dp*[1 0],T))/(2*dp);
% dy_2 = (idModels.func.fun_logdT([u1' u2'],p+dp*[0 1],T) - idModels.func.fun_logdT([u1' u2'],p-dp*[0 1],T))/(2*dp);
% subplot(2,1,1); plot(dy_1-dy(:,1)); subplot(2,1,2); plot(dy_2-dy(:,2));
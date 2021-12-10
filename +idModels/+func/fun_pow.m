function [F,dFdP] = fun_pow(u,p)
% Power function:
%
% F = mean(u,2)^p

u = mean(u,2);
F = u.^p;
dFdP = log(u).*F;
dFdP(u==0) = 0;
end

% x = 0:.001:1; dp = 1e-8; p = .75;
% [y,dy] = idModels.func.fun_pow(x',p);
% dy_ = (idModels.func.fun_pow(x',p+dp) - idModels.func.fun_pow(x',p-dp))/(2*dp);
% subplot(2,1,1); plot(x,y); subplot(2,1,2); plot(x,dy-dy_);
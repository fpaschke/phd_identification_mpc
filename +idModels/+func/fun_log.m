function [F,dFdP] = fun_log(u,p)
% Logarithmic function:
%
% F = ln(exp(p1)*mean(u,2) + 1)

u = mean(u,2);
F = log(exp(p(1))*u + 1);
dFdP = exp(p(1))*u./(exp(p(1))*u + 1);
end

% x = 0:.01:1; dp = 1e-8; p = -1e1;
% [y,dy] = idModels.func.fun_log(x',p);
% dy_ = (idModels.func.fun_log(x',p+dp) - idModels.func.fun_log(x',p-dp))/(2*dp);
% subplot(2,1,1); plot(x,y); subplot(2,1,2); plot(x,dy-dy_);
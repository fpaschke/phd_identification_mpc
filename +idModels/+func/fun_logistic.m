function [F,dFdP] = fun_logistic(u,p)
%Logistic function 
%
%                   1
% f(u) = --------------------------------
%           1 + exp( -p(2)*(u - p(1)) ) 
%
%u(1,:) -> the independent variable
%p(1) -> parameter to shift logistic function to left/right
%p(2) -> parameter to scale logistic fun 

x = u(:);
k = p(2);
x0 = p(1);
t = exp(-k*(x-x0));
F = 1./(1+t);
dFdP = [-k*t  (x-x0).*t]./(1+t).^2;
end

% idModels.test.testGrad(@idModels.func.fun_logistic,linspace(0,100,100)',[50 .1]);
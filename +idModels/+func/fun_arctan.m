function [F,dFdP] = fun_arctan(u,p)
%Arcus tangens 
%                   
% f(u,p) = arctan( k(x-x0) )/pi + 1/2
%
%u(1,:) -> the independent variable
%p(1) -> parameter to shift left/right
%p(2) -> parameter to scale

k = p(2);
x0 = p(1);
t = k*(u-x0);
F = atan(t)/pi +.5;
dFdP = [-k./(1+t.^2) (u-x0)./(1+t.^2)]./pi;
end

% x = 0:100; dp = 1e-8; p = [25 -1];
% [y,dy] = idModels.func.fun_arctan(x',p);
% dy_1 = (idModels.func.fun_arctan(x',p+dp*[1 0]) - idModels.func.fun_arctan(x',p-dp*[1 0]))/(2*dp);
% dy_2 = (idModels.func.fun_arctan(x',p+dp*[0 1]) - idModels.func.fun_arctan(x',p-dp*[0 1]))/(2*dp);
% subplot(2,1,1); plot(x,y); subplot(2,1,2); plot(x,dy-[dy_1 dy_2]);
function [F,dFdP] = fun_poly_normalized_HPeff(u,p,y)
X = u(:,1);
[F,dFdP] = idModels.func.fun_poly_normalized([y u(:,2)],p(2:end));
F = F.*X/p(1);
dFdP = [-F/p(1) dFdP.*X/p(1)];
end

% idModels.test.testGrad(@idModels.func.fun_poly_normalized_eff,[linspace(0,100,31)' (0:30)'],[100 10 3 6 1 1 1 1],(6:36)')
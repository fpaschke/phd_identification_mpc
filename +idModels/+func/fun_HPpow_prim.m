function [F,dFdP] = fun_HPpow_prim(u,p,y)

% u(:,1)    ... Ansteuerung
% u(:,2)    ... Tprim.sup
% y         ... Tsec,sup
P0 = p(1);  % Konstanter Leisungsverlust
a = p(2);   % Leistungsfaktor
[COP,dCOP] = idModels.func.fun_poly_normalized(y-u(:,2),p(3:6));
[Psec,dPsec] = idModels.func.fun_HPpow_sec(u(:,1:2),p(7:end),y);
F = ((COP-a)./COP).*Psec + (u(:,1)>0)*P0;
dFdP = [(u(:,1)>0) -Psec.*COP./COP.^2 a*Psec.*dCOP./(COP.^2) (COP-a)./COP.*dPsec];

% idModels.test.testGrad(@idModels.func.fun_HPpow_prim,[(0:100)' (0:.1:10)'],[.1 .9 5 35 1 .5 ones(1,8)],(30:.2:50)')
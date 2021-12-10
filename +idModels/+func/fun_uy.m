function [F,dFdP] = fun_uy(u,p,y)
% Outputs p(1)*y if u>0 and p(2)*y otherwise.
F = p(1)*(u>0).*y + p(2)*(u<=0).*y;
dFdP = [(u>0).*y (u<=0).*y];
end
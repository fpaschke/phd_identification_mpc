function [F,dFdP] = fun_y(u,p,y)
% p(1) == 1: Outputs y if u>p and 0 otherwise.
% p(1) == 0: Outputs y if u<=p and 0 otherwise.

if p(1) == 1
    F = (u>0).*y;
elseif p(1) == 0
    F = (u<=0).*y;
end
dFdP = NaN(length(F),1); % p1 is not meant for optimization
end
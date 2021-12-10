function [F,dFdP] = fun_Tank_Cond(u,p,y)
    %y ... T_upperlayer,T_thislayer,T_lowerlayer

N = size(y,1);
if size(y,2) == 2 %upper or lower layer
    F = y(:,2) - y(:,1);
else %mid layer
    F = y(:,1) + y(:,3) - 2*y(:,2); 
end
dFdP = NaN(N,0);
end
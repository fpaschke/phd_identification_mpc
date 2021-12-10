function [F,dFdP] = fun_Tank_Up(u,p,y)
    %u ... Vprim,Vsec
    %y ... T_thislayer,T_lowerlayer
    
dV = u(:,2) - u(:,1); %Vsec-Vprim
N = length(dV);
F = (dV>0).*dV.^p(1).*(y(:,2) - y(:,1));
dFdP = NaN(N,1);

end
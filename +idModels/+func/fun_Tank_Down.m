function [F,dFdP] = fun_Tank_Down(u,p,y)
    %u ... Vprim,Vsec
    %y ... T_upperlayer,T_thislayer
    
dV = u(:,1) - u(:,2); %Vprim-Vsec
N = length(dV);
F =  (dV>0).*dV.^p(1).*(y(:,1)-y(:,2));
dFdP = NaN(N,1);
end
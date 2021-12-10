function [F,dFdP] = fun_Tank(u,p,y)
    %u ... Vprim,Vsec,Tprim/Tsec
    %y ... Tupper,Tm,Tlower
    
dV = u(:,2) - u(:,1); %Vsec-Vprim > 0 -> unten nach oben 
switch p(1)
    case 1 %Upper layer
        %F = u(:,1).*(u(:,3) + p(2) - y(:,1)) + (dV>0).*dV.*(y(:,2) + p(3) - y(:,1));
        %dFdP = [NaN(length(dV),1) u(:,1) (dV>0).*dV];
        F = u(:,1).*(u(:,3) + p(2) - y(:,1)) + u(:,2).*(y(:,2) + p(3) - y(:,1));
        dFdP = [NaN(length(dV),1) u(:,1) u(:,2)];
    case 2 %Middle layer
        %F = (dV>0).*dV.*(y(:,3) + p(2) - y(:,2)) - (dV<0).*dV.*(y(:,1) + p(3)  -y(:,2));
        %dFdP = [NaN(length(dV),1) (dV>0).*dV  -(dV<0).*dV];
    	F = u(:,2).*(y(:,3) + p(2) - y(:,2)) + u(:,1).*(y(:,1) + p(3)  - y(:,2));
        dFdP = [NaN(length(dV),1) u(:,2) u(:,1)];
    case 3 %Lower layer
        %F = u(:,2).*(u(:,3)+ p(2) - y(:,2)) - (dV<0).*dV.*(y(:,1) + p(3) - y(:,2));
        %dFdP = [NaN(length(dV),1) u(:,2) (dV<0).*dV];
     	F = u(:,2).*(u(:,3)+ p(2) - y(:,2)) + u(:,1).*(y(:,1) + p(3) - y(:,2));
        dFdP = [NaN(length(dV),1) u(:,2) u(:,1)];
end

end
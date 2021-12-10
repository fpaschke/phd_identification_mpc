function [F,dFdP] = fun_VdT_switched(u,p,y)
% V multiplied by dT
% u(:,1) -> V_1
% u(:,2) -> V_2
% u(:,3) -> Tupper (mode 2)/Tsup (mode 1)  
% u(:,4) -> Tlower (mode 2)
% y(:,1) -> Tsp
% y(:,2) -> Tin (mode 1)/Tupper (mode 2)
% y(:,3) -> Tlower (mode 3)
% MODE 1:
%   (V_1 - V_2)*(Tin - Tsp)^p(1) if V_1 > V_2 else 0
% MODE 2:
%   |V_1 - V_2|*(Tupper - Tsp) if V1 > V2
%   |V_1 - V_2|*(Tlower - Tsp) if V2 < V1

V1 = u(:,1);
V2 = u(:,2);
V = V1 - V2;
Tsp = y(:,1);
ixV1grV2 = V>0;
if size(u,2)==4 || size(y,2) == 3
    mode = 2;
else 
    mode = 1;
end

if mode == 1
   	V(~ixV1grV2) = 0;
    if size(y,2) == 1
        Tsup = u(:,3);
    elseif size(y,2) == 2
        Tsup = y(:,2);
    else
        error('Unspecified yet!');
    end
else 
    V = abs(V);
 	if size(y,2) == 1
        Tsup = u(:,3).*ixV1grV2 + u(:,4).*~ixV1grV2;        
    elseif size(y,2) == 3
        Tsup = y(:,2).*ixV1grV2 + y(:,3).*~ixV1grV2;        
    else
        error('Unspecified yet!');
  	end
end
[f, df] = idModels.func.fun_dT(Tsup,p(1),Tsp);
F = V.*f;
dFdP = V.*df;
end
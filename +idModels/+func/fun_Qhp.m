function [F,dF] = fun_Qhp(u,p,y)
% p(end)=1 (sec): F = x*eps/(eps-1)*mprim*(TprimSup-TprimRet)
% p(end)=0 (prim): F = x*(eps-1)/(eps)*mprim*(TprimSup-TprimRet)
% u: Vprim, TprimSup, X     /       Vsec, TprimSup, X, TsecRet
% y: TsupSec, TretPrim      /       TsupSec, TretPrim 

[Eps,dEps] = idModels.func.fun_poly(y(:,1)-u(:,2),p(2:end-1)); 
if p(end) == 1 % sec
    [V, dV] = idModels.func.fun_VdT(u(:,[1 2]),p(1),y(:,2));
    Eta = Eps./(Eps-1);
    dEta = -dEps./(Eps-1).^2;
elseif p(end) == 0 % prim
    [V, dV] = idModels.func.fun_VdT([u(:,1) y(:,1)],p(1),u(:,4));
    Eta = (Eps-1)./Eps;
    dEta = dEps./Eps.^2;
else
    error('Last Parameter needs to be either 0 or 1');
end
F = V.*u(:,3).*Eta;
dF = u(:,3).*[dV.*Eta V.*dEta NaN(length(F),1)];
end

% idModels.test.testGrad(@idModels.func.fun_Qhp,10*rand(100,4),[1 10 1 1 0],zeros(100,2)) 

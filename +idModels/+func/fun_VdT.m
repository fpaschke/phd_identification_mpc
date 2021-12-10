function [F,dFdP,dFdU,dFdY] = fun_VdT(u,p,y)
% V multiplied by dT
% sum(V,2) * (Tsup - Troom)^p(1) 

if nargin < 3
    nv = size(u,2) - 2;
    Tsup = u(:,end-1);
    Tr = u(:,end);
else
    nv = size(u,2) - 1;
    Tsup = u(:,end);
    Tr = y(:,1);
end
if nv == 0 % Volumeflow/Massflowrate is determined by output 2!
    V = y(:,2);
else
    V = sum(u(:,1:nv),2);
end

switch nargout
    case 1
      	dT = idModels.func.fun_dT(Tsup,p(1:end),Tr);
    case 2
        [dT, dT_dP] = idModels.func.fun_dT(Tsup,p(1:end),Tr);       
    otherwise
       [dT, dT_dP, dTdU, dTdY] = idModels.func.fun_dT(Tsup,p(1:end),Tr);
end
F = V.*dT;
if nargout > 1; dFdP = V.*dT_dP; end
if nargout > 2
    N = size(u,1);
    dFdU = [ones(N,nv).*dT V.*dTdU];
    dFdY = V.*dTdY;
    if nv == 0
        dFdY = [dFdY ones(N,1).*dT];
    end
end
end

% idModels.test.testGrad(@idModels.func.fun_VdT,[(0:.1:2)' (10:1:30)'],[1.3 2 3],5*rand(21,1)+20)
% idModels.test.testGrad(@idModels.func.fun_VdT,[2*rand(20,1) 5*rand(20,1)+20],[1],5*rand(20,1)+20,'ix_out',3,'ix_p',1)
% idModels.test.testGrad(@idModels.func.fun_VdT,[2*rand(20,1) 5*rand(20,1)+20],[1],5*rand(20,1)+20,'ix_out',4,'ix_p',3)
function [F,dFdP,dFdU,dFdY] = fun_dT(u,p,y)
% Function:
%   y = sgn(dT)*abs(dT)^p(1) where 
% nargin == 3
%   dT = p(3)*u-p(2) - y
% nargin == 2
%   dT = p(3)*u(:,1)-p(2) - u(:,2)

if nargin < 3
    Ts = u(:,1);
    Tr = u(:,2);
else
    if isempty(u)
        Ts = y(:,1);
        Tr = y(:,2);
    else
        Ts = u(:,1);
        Tr = y;
    end
end

switch length(p)
    case 1
        dT = Ts - Tr;
    case 2
        dT = Ts - p(2) - Tr;
    case 3
        dT = p(3)*Ts - p(2) - Tr;
end

F = sign(dT).*abs(dT).^p(1);    % negative zahlen hoch ungerader exponent sind komplex!  
F(dT==0) = 0;                   % falls p<0
if nargout > 1
    N = size(Ts,1);
    dFdP(:,1) = log(abs(dT)).*F; % f'(a^x)=ln(a)*a^x
    dFdP(dT==0,1) = 0;
    if length(p)>1
        ind = zeros(N,1);
        ind(dT>0) = 1;
        ind(dT<0) = -1;
        x = p(1)*sign(dT).*ind.*abs(dT).^(p(1)-1);
        dFdP(:,2) = -x;
    end
   	if length(p) > 2
        dFdP(:,3) = Ts.*x;
    end
end
if nargout > 2 
    if mod(p(1),1)~=0 && nargin == 3 % Not implemented yet
        dFdU = NaN(size(F,1),1);
        dFdY = NaN(size(F,1),1);
    else
        if length(p) == 3; dTdU = p(3); else; dTdU = 1; end
        dFdU = p(1)*dTdU*dT.^(p(1)-1);
        dFdY = -p(1)*dT.^(p(1)-1);
    end
end
end
% x = -10:0.1:10; dp = 1e-8; p = [1.3 0 1]; T = zeros(length(x),1);
% [y,dy] = idModels.func.fun_dT(x',p,T);
% ixp = 3;
% pm = p; pp = p; pm(ixp) = pm(ixp) - dp; pp(ixp) = pp(ixp) + dp; 
% dy_ = (idModels.func.fun_dT(x',pp,T) - idModels.func.fun_dT(x',pm,T))/(2*dp);
% figure; subplot(2,1,1); plot(x,y); subplot(2,1,2); plot(x,dy(:,ixp)-dy_);
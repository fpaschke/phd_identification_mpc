function [F,dFdP] = fun_softplus(u,p)
% Softplus function:
%   y =   p(1)*ln(1 + exp[p(2) * ( u - p(3) ) ] )
% u -> the independent variable
% p Parametervector

u = u(:);
Npcs = length(p)/3;
F = 0; dFdP = []; np = 0;
for j = 1:Npcs
    k = p(np + 1);
	ku = p(np + 2);
    u0 = p(np + 3);
    
    t_0 = exp(ku*(u - u0));
    t_1 = 1 + t_0;
    F = F + k*log(t_1);
    if nargout > 1
        dFdP = [dFdP log(t_1) k*(u - u0).*t_0./t_1 -k*ku*t_0./t_1 ];
    end
    np = np + 3;
end
end

% T = -20:40; figure; plot(T,fun_softplus(T,[1 -1 15]))
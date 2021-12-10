function [F,dFdP] = fun_hcv(u,p)
% y = p1 (p2 - Tout)^p3 + p2 if Tout < p(2) otherwise y = p(2)
x = p(2)-u;
ix = x<=0;
x(ix) = 0;
F = p(1)*x.^p(3) + p(2);
dFdP = [x.^p(3) ...
        p(1)*p(3)*x.^(p(3)-1)+1 ...
        p(1)*log(x).*x.^p(3)]; % f'(a^x)=ln(a)*a^x
dFdP(ix,:) = [0 1 0].*ones(sum(ix),1);
end

% close all; idModels.test.testGrad(@idModels.func.fun_hcv,(-20:1:25)',[2 22 1/1.3],'delta',1e-8)
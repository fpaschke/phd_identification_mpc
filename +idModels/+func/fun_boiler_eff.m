function [F,dFdP] = fun_boiler_eff(u,p,y)
%Scalable Logistic function with offset
%
% IF u has two columns Tflow = u(:,2) or if y is supplied Tflow = y:
%
%                                           p(5)
% f(u) = u(:,1)/p(2)*p(1)(1 -	-------------------------------- )
%                               1 + exp( -p(4)*( Tflow - p(3) ) ) 
%
% ELSE
%
% f(u) = u(:,1)/p(2)*p(1)

u1 = u(:,1)/p(2); % Normierung des Ansteuersignals auf [0...1]
if size(u,2) == 1 && nargin <=2
    F = p(1)*u1;
    dFdP = [u1 -u1*p(1)/p(2)];
elseif size(u,2) == 2 || nargin > 2
    if nargin > 2
        Tflow = y;
    else
        Tflow = u(:,2);
    end
    [F,dFdP] = idModels.func.fun_logistic(Tflow,p(3:4));
    x = (1 - F*p(5));
    dFdP = [u1.*x -u1.*x*p(1)/p(2) -p(1)*u1.*[p(5)*dFdP F]];
    F = p(1)*u1.*x;
end
end

% idModels.test.testGrad(@idModels.func.fun_boiler_eff,[100*ones(1,100); linspace(1,100,100)]',[285 100 50 1 .2]);
function [f,df] = fun_radlogistic(u,p)
%calculates effective radiation through the windows weighted with logistic function 
%
%   f(u) = fun_logistic(u,p) * fun_rad(u,p)
%
%u(1,:) -> global rad
%u(2,:) -> azimuth
%u(3,:) -> sun height
%p(1) -> rotation angle from north (west->270) (fixed parameter!)
%p(2) -> horizontal orientation (perpendicular to surface of earth ->90) (fixed parameter!)
%p(3) -> parameter to shift logistic function to left/right
%p(4) -> factor for logistic fun 

[W,dW] = idModels.func.fun_rad(u,p(1:2));
[L,dL] = idModels.func.fun_logistic(u(:,3),p(3:4));
f = W.*L;
df = [dW.*L dL.*W];

% if any(isnan(f)) || any(~isreal(f)) || any(isnan(dL(:))) || any(~isreal(dL(:)))
%    stop = 1; 
% end

end
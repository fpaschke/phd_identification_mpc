function [f,df] = fun_radarctan(u,p)
%calculates effective radiation through the windows weighted with arctan function 
%
%   f(u) = fun_arctan(u,p) * fun_rad(u,p)
%
%u(1,:) -> global rad
%u(2,:) -> azimuth
%u(3,:) -> sun height
%p(1) -> rotation angle from north (west->270) (fixed parameter!)
%p(2) -> horizontal orientation (perpendicular to surface of earth ->90) (fixed parameter!)
%p(3) -> parameter to shift arctan function to left/right
%p(4) -> steepnessfactor of arctan fun 

[W,dW] = idModels.func.fun_rad(u,p(1:2));
[L,dL] = idModels.func.fun_arctan(u(:,3),p(3:4));
f = W.*L;
df = [dW.*L dL.*W];

if any(isnan(f)) || any(~isreal(f)) || any(isnan(dL(:))) || any(~isreal(dL(:)))
   stop = 1; 
end

end
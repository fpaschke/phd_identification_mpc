function [f,df] = fun_rad_dir(u,p)
%calculates effective radiation through the windows 
%u(1,:) -> direct normal radiation 
%u(2,:) -> azimuth
%u(3,:) -> sun height
%u(4,:) -> shading position (if size(u,2)>3) 0 ... opened, 1 ... closed
%p(1) -> rotation angle from north (west->270) 
%p(2) -> horizontal orientation (perpendicular to surface of earth ->90)

l_h = 180-u(:,3)-p(2); % angle between sun height and plane 
l_az = 90-p(1)+u(:,2); % angle between sun azimuth and plane [-360 ... 360]
f = u(:,1).*sin(l_h*pi/180).*sin(l_az*pi/180);
f(u(:,3)<0 | sin(l_az*pi/180)<0)  = 0; % sun below horizon or behind window
if size(u,2)>3
    f = (1-u(:,4)).*f;
end
if nargout > 1 
    df = NaN(length(f),2); %Keine Lust Gradienten zu berechnen... Bisher auch noch nicht nötig, da keine Optimierungsparameter
end
end
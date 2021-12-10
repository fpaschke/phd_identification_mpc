function [f,df] = fun_rad_diff(u,p)
%calculates radiation through the windows 
%u(1,:) -> diff radiation
%u(2,:) -> shading position (if size(u,2)>3) 0 ... opened, 1 ... closed

f = u(:,1);
if size(u,2)>1
    f = (1-u(:,2)).*f;
end
if nargout > 1 
    df = NaN(length(f),0);
end
end
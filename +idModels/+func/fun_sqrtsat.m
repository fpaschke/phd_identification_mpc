function F = fun_sqrtsat(u,p)
% Squarroot Function with or without Saturation
%   y =   p1* (sqrt(x).*(x<=p2) + sqrt(p2).*(x>p2)) , where x = u1 + u2 + ... up
% u -> the independent variable
% p Parametervector

x = sum(u,2);
if length(p) > 1
    F  =   p(1)* (sqrt(x).*(x<=p(2)) + sqrt(p(2)).*(x>p(2)));
else
    F = p(1)*sqrt(x);
end
end
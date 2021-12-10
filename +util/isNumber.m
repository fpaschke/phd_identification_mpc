function y = isNumber(x)
%ISNUMBER Check weather x is finite ~NaN number
%
% y = isNumber(x)
% x [double array] 
% y [int array]
%
% o = isNumber([1 Inf NaN -Inf])

if iscell(x)
    y = false;
else
    y = isnumeric(x) & ~isnan(x) & ~isinf(x);
end


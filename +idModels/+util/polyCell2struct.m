function Ps = polyCell2struct(Pc)
%POLYCELL2STRUCT Convert Cell Polynomial to struct.
%
% Ps = polyCell2struct(Pc)
% Ps [n x m struct]: Struct with field den and num.
% Pc [n x m cell]: cell of 2 x 1 cells where first element defines numerator and second element defines denominator.

assert(iscell(Pc),'Input needs to be a nr x nc cell of 2 x 1 cells of double vectors!');
if isnumeric(Pc{1}) && isvector(Pc{1})
    Pc = {Pc};
end
    
[nr,nc] = size(Pc); 
for r = 1:nr
    for c = 1:nc
        Ps(r,c).num = Pc{r,c}{1}(:)';
        Ps(r,c).den = Pc{r,c}{2}(:)';
    end
end
end


function isDiag = isDiagPoly(P)
%ISDIAGPOLY Checks if polynomial matrix is an diagonal Poly Matrix.
%All fields (val,free,std,factorized) of corresponding poly need to be equal!
%
% P [n x m struct array]: Struct representation of polynomial matrix
% isDiag [logical]: If true the polynomial matrix is diag otherwise not

isDiag = true;
[Nr,Nc] = size(P);
if Nr ~= Nc
    isDiag = false;
    return;
end
for nr = 1:Nr
    for nc = 1:nr-1
        if nr ~= nc
           if ~isequaln(P(nr,nc),P(nc,nr))
               isDiag = false;
               break;
           end
        end
    end
end
end


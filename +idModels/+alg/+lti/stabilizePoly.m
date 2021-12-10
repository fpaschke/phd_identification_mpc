function P = stabilizePoly(P,fac,tr)
% Stabilizes Polynomial P by deviding unstable roots by their manitude*fac
%
% Ps = stabilizePoly(P,fac)
% Ps [double vector]: Stabilized Polynomial coefficients 
% P [double vector]: Polynomial coefficients to be stabilized 
% fac [opt. double scalar]: Factor by which unstable modes will be devided
% tr [opt. double scalar]: threshold for which poles are regarded as unstable

if any(isnan(P.val))
    return;
end
if nargin <= 1 || isempty(fac)
    fac = 1.1;
end
if nargin <= 2 || isempty(tr)
    tr = 1;
end
assert(isscalar(fac) && fac >= 1,'"fac" needs to be a scalar >1');
if P.factorized
    r_c = P.val(2:end);
else
    r_c = roots(P.val);
end
mag_rc = abs(r_c);
ix_unstab = mag_rc>=tr;
if any(ix_unstab)    
    r_c(ix_unstab) = r_c(ix_unstab)./(mag_rc(ix_unstab)*fac);
    if ~P.factorized
        P.val = poly(r_c);
    else
        P.val(2:end) = r_c;
    end
end 
end
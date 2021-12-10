function P0 = factorizePoly(P)
%Computes zeros p0i from polynomial coefficients pi of P(q):
% P(q) = p0 + p1 q^-1 + ... + pn q^-n = K*(1 - p10 q-^1) * ... * (1 - pn0 q-^1)
% A zero at Inf corresponds to a delay.
%
% P0 = factorizePoly(P,imagTol)
% P [n+1 x 1 double]:  Vector of polynomial coefficients [p0 p1 ... pn]
% P0 [n+1 x 1 double]: Vector representing Gain and poles of P [K p01 ... p0n]. A zero at Inf corresponds to a delay.

if isempty(P)
    P0 = P;
elseif all(P==0)
    P0 = [0 Inf(1,length(P)-1)];
else
    nk = find(P~=0,1,'first')-1;
	K = P(nk+1);
    if all(~isnan(P(nk+1:end)))
        pr = roots(P(nk+1:end)/K);
    else
        pr = NaN(length(P(nk+1:end))-1,1);
    end
    P0 = [K Inf(1,nk) pr'];   
end
end
function [P,dP_dp0] = defactorizePoly(P0,imagTol)
%Computes polynomial coefficients pi from zeros p0i of P(q):
% P(q) = p0 + p1 q^-1 + ... + pn q^-n = K*(1 - p10 q-^1) * ... * (1 - pn0 q-^1)
%and the derivatives of pi wrt to the zeros.
%
% [P,dP_dp0] = defactorizePoly(A,imagTol)
% P0 [n+1 x 1 double]:  Vector representing Gain and zeros of P [Kp p01 ... p0n]. Zero at Inf corresponds to delay.
% P [n+1 x 1 double]: Vector representing coeffitients pi [p0 ... pn]
% dP_da [n+1 x n+1 matrix]: Vector representing derivatives of pi wrt Gain K and poles p0i [dp0/dK ... dpn/dK; ...; dp0/dp0n .... dpn/dp0n]

if nargin == 1
    imagTol = 1e-12;
end

if isempty(P0)
    P = [];
    dP_dp0 = [];
else
    n = length(P0) - 1;
    ix_delay = P0(2:end)==Inf;
    rts = P0(find(~ix_delay)+1);
    if all(~isnan(rts))
        R = [zeros(1,sum(ix_delay)) poly(rts)];
        P = R*P0(1);
        ix = abs(imag(P))<=imagTol;
        P(ix) = real(P(ix));
        R(ix) = real(R(ix));
    else
        P = [zeros(1,sum(ix_delay)) P0(1) NaN(1,length(rts))];
        R = P/P0(1);
    end

    if nargout > 1
        dP_dp0 = zeros(n+1,n+1);
        dP_dp0(1,:) = R; 
        for np = 2:n+1
            dP_dp0(np,2:end) = -P0(1)*deconv(R,[1 -P0(np)]);
        end
    end
end
end
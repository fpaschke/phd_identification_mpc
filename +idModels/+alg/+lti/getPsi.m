function psi = getPsi(P,N)
% GETPSI Calculates Matrix only PSI such that Y = PSI*y0 + THETA*U, where 
% PSI = [C, CA, ... CA^Hp]'; 
% THETA = [D 0 ...  0; CB D 0 ...0; CAB CB D 0 ... 0; ... ; CA^Hp-1B ...]'
% of the Statespace description of P(q^-1)*y(k) = u(k) 
% which can be used to calculate initial conditions of the filter y = 1/P*u
% 
% psi = getPsi(P,N)
% PSI [N x length(P)-1 double]: Extended observability Matrix.
% P [N x length(P)-1 double or matlab tf object]: monic polynomial or matlab tf object
% PSI [N x length(P)-1 double]: Extended observability Matrix.

assert(isscalar(N) && mod(N,1)==0 && N>0,'"N" should be a pos. scalar!');
if isvector(P) && isnumeric(P)
    if iscell(P)
        P = P{1,1};
    end
    assert(isvector(P) && isnumeric(P),'"P" should be a numric vector!');
    assert(P(1) == 1,'"P" should be a monic polynomial!');

    n = length(P) - 1;
    A = diag(ones(n-1,1),-1); A(1,:) = -P(2:end); C = [1; zeros(n-1,1)]';
    if nargout > 1
       G = tf(1,P,1,'Variable','z^-1');
    end
else
    Gss = minreal(ss(P),[],false);
    A = Gss.A; C = Gss.C;
    n = length(A);
end
ny = size(C,1);
psi = zeros(N*ny,n); psi(1:ny,:) = C*A;
for k = 2:N
    r = (k-1)*ny+1:k*ny;
    psi(r,:) = psi(r-ny,:)*A;
end
end
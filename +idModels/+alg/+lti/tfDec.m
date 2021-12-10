function [F, Bd] = tfDec(B,A,k)
%TFDEC Decompose transfer function G(q^-1) = B(q^-1)/A(q^-1) to 
% F(q^-1) + Bd(q^-1)/A(q^-1), where F = h(0) + h(1)q^-1 + ... + h(k)q^-k 
% contains the first k+1 coefficients of the series expansion of G.
%
% [F, Bd, Ad] = tfDec(B,A,l)
% B [double vector]: Coeff. of numerator polynomial
% A [double vector]: Coeff. of denom polynomial
% k [scalar int]: number of coeff. of F
% F [1 x k+1 double vector]: coefficients [h(0) h(1) ... h(k-1)] of F
% Bd [double vector]: coefficients of "new" numerator polynomial 

assert(isvector(B) && isvector(A) && isnumeric(A) && isnumeric(B) && isreal(A) && isreal(B));
assert(isscalar(k) && k>=0 && mod(k,1)==0);
A = A(:)'; B = B(:)';
k = k + 1;

na = length(A)-1; nb = length(B)-1;
n = max(na,nb);
ng = max(na-1,nb-k);
C_ = [B zeros(1,n-nb)]; 
A_ = [A zeros(1,n-na)]; 
[F,Bd] = deconv([C_ zeros(1,k-1)],A_);
end

% a = conv([1 -.9],[1 -.8]); b = [1 .5 .3 .4 .5];
% [F, Bd] = idModels.alg.lti.tfDec(b,a,0)
% x = randn(1e3,1);
% sum(abs(filter(b,a,x) - filter(F,1,x) - filter(Bd,a,x)))
 


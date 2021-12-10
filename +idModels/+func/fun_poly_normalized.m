function [F,dFdP] = fun_poly_normalized(u,p)
% Polynomial
%   y = p0 + p11 (u_1 - u10) + ... + p_1n (u_1 - u10)^n + ... + p_nu_n (u_n - u_n0) u^n
% n = (length(p)-nu-1)/nu
% p = (F0 u_10 ... u_n0 p11 ... p_1n ... p_nu_1 ... p_nu_n)

p = p(:);               % Make sure that p is an rowvector
[N,n_u] = size(u);      % Samples/Number of Ins
F0 = p(1);              % Intercept
u0 = p(2:n_u+1)';       % Normalization
p = p(n_u+2:end);       % Coefficients
un = u - u0;            % Normalized Inputs
np_u = length(p)/n_u;   % Degree of Polynomials
assert(mod(np_u,1)==0); 

dFdP = [ones(N,1)  NaN(N,length(p)+n_u)]; % Init Gradient
ixp = 1;
for i = 1:n_u
    dFdCoeff = bsxfun(@power,un(:,i),1:np_u);
    p_nu = p(ixp:ixp+np_u-1)';
    dFdP(:,i+1) = -p_nu(1)-sum(p_nu(2:end).*(2:np_u).*dFdCoeff(:,1:end-1),2);  % Gradient wrt normalization Parameters
    dFdP(:,1+n_u+ixp:n_u+ixp+np_u) = dFdCoeff;
    ixp = ixp + np_u;
end
F = F0 + dFdP(:,n_u+2:end)*p;
end

% idModels.test.testGrad(@idModels.func.fun_poly_normalized,[(0:30)' (6:36)'],[10 3 6 1 1 1 1])
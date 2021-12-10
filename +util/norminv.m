function x = norminv(p,mu,sigma)
%NORMINV This is an workaround if matlabs norminv function  is not available.
% 
% x = norminv(p,mu,sigma)
% p [double vector]: Probability values 0...1
% mu [opt double scalar]: Mean. Def. 0
% sigma [opt double scalar]: Standard dev. Def. 0

% TODO: Call Matlabs norminv if stats toolbox is available

if nargin<2 || isempty(mu)
    mu = 0;
end
if nargin<3 || isempty(sigma)
    sigma = 1;
end

assert(isvector(mu) && isvector(sigma) && length(mu) == length(sigma) && (isscalar(p) || (length(p) == length(mu)) ));
x = sigma.*(-sqrt(2)*erfcinv(2*p))+mu;
end


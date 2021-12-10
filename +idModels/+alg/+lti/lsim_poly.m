function y = lsim_poly(A,B,C,u,e,y0,u0,e0)
%LSIM_POLY Simulate discrete time mimo system 
% y(k) + A1 y(k-1) + ... + Ana y(k-na) =  B0 u(k) + ... + Bnb u(k-nb) + e(k) + C1 e(k-1) + ... + Cnc e(k-nc)
%
% y = lsim_poly(A,B,C,u,e,y0,u0,e0)
% 
% EXAMPLE:
% A = {[-.95 -.1; 0 -.9]}; ny = size(A{1},2); nu = 2; B = {zeros(ny,nu) eye(ny,nu)};  C = {eye(ny)};
% N = 100; u = rand(N,nu); e = .1*randn(N,ny);
% y = idModels.alg.lti.lsim_poly(A,B,C,u,e);

% Init
if ~isempty(u); Ns = size(u,1); else; Ns = size(e,1); end 
ny = size(A{1},1);
nu = size(u,2);
if nargin <= 2
    C = [];
end
if ismatrix(A) && ~iscell(A)
	A = {A};
end
if ismatrix(B) && ~iscell(B)
	B = {B};
end
if isempty(C)
    C = cell(1,0);
end
if ismatrix(C) && ~iscell(C)
	C = {C};
end
if nargin <= 4 || isempty(e) % if e not supplied then create empty e
    e = zeros(Ns,0);  
end
if isempty(u) % if u not supplied then create empty u
    u = zeros(Ns,0);  
end
na = length(A);
nb = length(B)-1;
nc = length(C);
dmax = max([na nb nc]);

% build Theta
if isempty(e)
    Theta = [-cell2mat(A) cell2mat(B)]';
else
    Theta = [-cell2mat(A) cell2mat(B) eye(ny) cell2mat(C)]';
end

% Malloc
if nargin<=5 || isempty(y0)
    y0 = zeros(na,ny);
end
y = [nan(dmax-na,ny); y0; nan(Ns,ny)];

if nargin<=6 || isempty(u0)
    u0 = zeros(nb,nu);
end
u = [nan(dmax-nb,nu); u0; u];

if nargin<=7 || isempty(e0)
    e0 = zeros(nc,size(e,2));
end
e = [nan(dmax-nc,size(e,2)); e0; e];

% Sim
for k = dmax+1:Ns+dmax
    yy = y(k-1:-1:k-na,:)'; uu = u(k:-1:k-nb,:)'; ee = e(k:-1:k-nc,:)';
    phiT = [yy(:); uu(:); ee(:)]';
    y(k,:) = phiT*Theta;
end
y = y(dmax+1:end,:);
end
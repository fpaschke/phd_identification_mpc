function [F,dFdP] = fun_HPpow_sec(u,p,y)
% Polynomial

u0 = p(3:4);        % Normalization
[N,nu] = size(u);  	% Samples/Number of Ins
if nargin < 3     	
    % u = [TsupSec TsupPrim]
    if nu == 2
        un = u - u0;  	% Normalized Inputs
        X = ones(N,1);
    else % u = [X TretSec TsupPrim]
     	un = u(:,2:3) - u0;  	% Normalized Inputs
        X = u(:,1)/p(1);
    end
else
  	% u = [x TsupPrim]
    un = [y u(:,2)] - u0;  	% Normalized Inputs
    X = u(:,1)/p(1);
end
    F = X.*(p(2) + p(5)*un(:,1) + p(6)*un(:,2) + ...
             + p(7)*un(:,1).*un(:,2) + p(8)*un(:,1).^2 + p(9)*un(:,2).^2);
    dFdP = X.*[NaN(N,1) ones(N,1) NaN(N,2) un(:,1:2) un(:,1).*un(:,2) un(:,1:2).^2];
end

% idModels.test.testGrad(@idModels.func.fun_HPpow_sec,[(0:30)' (6:36)'],[100 10 3 6 1 1 1 1 1])
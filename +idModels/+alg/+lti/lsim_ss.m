function [y,x] = lsim_ss(A,B,C,D,u,x0)
%Simulate discrete linear statespace system.

Ns = size(u,1);
n = size(A,1);
ny = size(C,1);
x = NaN(n,Ns+1);
y = NaN(ny,Ns);

x(:,1) = x0;
for k = 1:Ns
    y(:,k) = C*x(:,k) + D*u(k,:)';
    x(:,k+1) = A*x(:,k) + B*u(k,:)';
end
y = y'; x = x(:,1:end-1)';
end


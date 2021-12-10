function [f,df_deta,df_du,df_dy] = f_2(u,eta,y)
% Berechnet Funktionswert f_2 und deren Ableitungen 
    f = u.^2.*y + eta*u.^3;           % f(u,y,eta) 
    df_deta = u.^3;                   % df_deta(u,y,eta)  
    df_du = 2*y.*u+3*eta.*u.^2;       % df_du(u,y,eta)
    df_dy = u.^2;                     % df_dy(u,y,eta)
end
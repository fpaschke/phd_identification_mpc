function [f,df_deta,df_du,df_dy] = f_1(u,eta,y)
% Berechnet Funktionswert f_1 und deren Ableitungen 
    f = sum(u,2).*(eta - y);          % f(u,y,eta)
    df_deta = sum(u,2);               % df_deta(u,y,eta)
    df_du = ones(1,2).*(eta - y);     % df_du(u,y,eta)
    df_dy = -sum(u,2);                % df_dy(u,y,eta)
end


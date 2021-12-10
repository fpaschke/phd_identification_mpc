clc; close all; clear;
addpath('..\..\..')

%% SETUP
% MA Poly
c_nC = -.99;                    % Parameter                   
nC = unique(floor(logspace(1,3,30))); % Orders of C Polynomial

% Data
rng(0);                         % For reproducability
N_s = 10e3;                     % Number of samples
N_exp = 1e2;                    % Number of experiments
e = randn(N_s,N_exp);           % Noise input to MA Model

%% Test
k = 1;
for nc = nC
    C = [1 zeros(1,nc-1) c_nC];
    w = filter(C,1,e); %MA Filtering
    M = idModels.NsfPolyModel(0,[],[],nc); %Model Setup
    M.C.val = C;
    for ne = 1:N_exp
        % BC
        tic;
        eps0 = M.icBackcast(w(:,ne));
        t_bc(k,ne) = toc;
        eps_bc = filter(1,C,w(:,ne),eps0);
        mse_bc(k,ne) = mean(eps_bc.^2);
        % LS
        tic;
        PSI = idModels.alg.lti.getPsi(C,N_s);
     	R = filter(1,C,w(:,ne));
      	eps0 = PSI\(-R); 
        t_ls(k,ne) = toc;
        eps_ls = PSI*eps0 + R;
        mse_ls(k,ne) = mean(eps_ls.^2);
    end 
    k = k + 1
end

%% Plot
figure('color','w','Position',[200 200 400 300]);
loglog(nC(1:size(t_bc,1)),[mean(t_ls,2) mean(t_bc,2)]); grid on; legend({'MkQ' 'Backforecasting'},'Location','NorthWest');
util.formatFigure(14,'Polynomordnung $n_L$','Rechenzeit in s');
util.saveTightFigure(gcf,'D:\Diss\Bilder\Anhang\IcEstTime.pdf')

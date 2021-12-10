clc; close all; clear;
addpath('..\..\..')

%% SETUP
% MA Poly
C = [1 zeros(1,0) -.99];            % C polynomial       
z0 = 10*ones(1,length(C));             % Act IC

% Data
rng(0);                             % For reproducability
N_exp = 100;                        % Number of experiments
%N_s = round(logspace(1,3,30));      % Number of samples
N_s = 20:20:1000;          
e = randn(max(N_s),N_exp);          % Noise input to MA Model

%% Test

for k=1:length(N_s)
    w = filter(C,1,e(1:N_s(k),:),1); %MA Filtering
    M = idModels.NsfPolyModel(0,[],[],length(C)-1); %Model Setup
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
        PSI = idModels.alg.lti.getPsi(C,N_s(k));
     	R = filter(1,C,w(:,ne));
      	eps0 = PSI\(-R); 
        t_ls(k,ne) = toc;
        eps_ls = PSI*eps0 + R;
        mse_ls(k,ne) = mean(eps_ls.^2);
        % Zero init
     	eps_zero = filter(1,C,w(:,ne));
        mse_zero(k,ne) = mean(eps_zero.^2);
    end 
    k = k + 1
end

%% Plot
figure('color','w','Position',[200 200 400 300]);
plot(N_s,[mean(mse_ls,2) mean(mse_bc,2) mean(mse_zero,2)]);  grid on; legend({'MkQ' 'Backforecasting' 'Anfangswert 0'},'Location','NorthEast');
set(gca,'XTick',0:200:1000)
util.formatFigure(14,'Datenpunktanzahl $N$','MSE($\varepsilon [t]$)'); 
util.saveTightFigure(gcf,'D:\Diss\Bilder\Anhang\IcEstPrec.pdf')

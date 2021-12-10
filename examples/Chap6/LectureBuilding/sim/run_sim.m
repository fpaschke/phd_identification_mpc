% Simulate temperature controlled building with 4 lecture halls.
% The simulation starts on P.Tstart at 00:00 which is assumed to be a Saturday
% A thorough description of the model can be found in the package
% LecBuilding.mo. 

clear; clc;   
addpath('..\..\..\..');
addpath('..\model\')
addpath('.\resource\')
addpath('..\..\SingleRoom\sim\subfun\ctrl\')
addpath('..\..\SingleRoom\sim\subfun\sig\')
addpath(genpath('.\subfun\'))
import idModels.*;

%% Create Parameters and Input Signals for Simulation
P = getParameters('T_s',5,'T_s_mpc',60);

%% Simulate PI Control Loop for different Heatup Times
P.CtrlFun = 'ctrl_pi';          % Name of controller function to be called 
Tpre = 0.5;
for k = 1:length(Tpre)
    P.Ctrl.PI.T_pre(:) = Tpre(k);
    SimOut_PI{k} = simulate('LecBuil',P);
end

%% Simulate MPC Control Loop with different settings
P.CtrlFun = 'ctrl_mpc';         % Name of controller function to be called 
P.N_w = 52;
P.Ctrl.MPC.MpcActivationTime = '2021-1-1';
P.Ctrl.MPC.UseWeatherPredictionModel = false;
P.Ctrl.MPC.UseMimoModel = false;
SimOut_MPC{1} = simulate('LecBuil',P);
P.Ctrl.MPC.UseMimoModel = true;
SimOut_MPC{2} = simulate('LecBuil',P);
P.Ctrl.MPC.UseWeatherPredictionModel = true;
P.Ctrl.MPC.UseMimoModel = false;
SimOut_MPC{3} = simulate('LecBuil',P);
P.Ctrl.MPC.UseWeatherPredictionModel = true;
P.Ctrl.MPC.UseMimoModel = false;
P.Ctrl.MPC.UseSeasonalModel = true;
P.Ctrl.MPC.UseSparse = true;
SimOut_MPC{4} = simulate('LecBuil',P);

%% Plot 
ix_pi = 1;
ix_mpc = 1;

% Plot Trajectories
figure('Color','w','Position',[100 100 900 1100]);
Nz = size(SimOut_PI{ix_pi}.m_sup.Data,2);
lcol = lines(2);
for nz = 1:Nz
    s(nz) = subplot(Nz+1,1,nz);
    plot(getSubSignal(SimOut_PI{ix_pi}.T_air,nz),'Color',lcol(2,:),'LineStyle','--','DisplayName','Luft Temp. (PI) [C]'); hold on;
    plot(getSubSignal(SimOut_MPC{ix_mpc}.T_air,nz),'Color',lcol(2,:),'LineStyle','-','DisplayName','Luft Temp. (MPC) [C]'); hold on;
    plot(getSubSignal(SimOut_PI{ix_pi}.In.T_ref_lower,nz),'--k','DisplayName','Unterer Grenzwert [C]'); hold on;
    plot(getSubSignal(SimOut_PI{ix_pi}.In.T_ref_upper,nz),'--k','DisplayName','Oberer Grenzwert [C]'); hold on;
    ylabel(['$\vartheta_' num2str(nz) ' [^\circ \mathrm{C}$]']); ylim([19 31]); yyaxis right; 
    plot(getSubSignal(SimOut_PI{ix_pi}.m_sup,nz),'Color',lcol(1,:),'LineStyle','--','DisplayName','Massenstrom (PI) [kg/s]'); hold on; 
    plot(getSubSignal(SimOut_MPC{ix_mpc}.m_sup,nz),'Color',lcol(1,:),'LineStyle','-','DisplayName','Massenstrom (PI) [kg/s]'); grid on;
    if nz < 3; ylim([-.5 6.5]); else; ylim([-.5 4.5]); end
    s(nz).YColor = lcol(1,:); ylabel(['$\dot{m}_' num2str(nz) '$ [kg/s]']); title(''); yyaxis left; s(nz).XTickLabel = {};
    yyaxis left; s(nz).YColor = lcol(2,:); 
end
s(Nz+1) = subplot(Nz+1,1,Nz+1);
plot(SimOut_PI{ix_pi}.T_sup,'Color',lcol(2,:),'LineStyle','--','DisplayName','Zulufttemp. (PI) [C]'); hold on;
plot(SimOut_MPC{ix_pi}.T_sup,'Color',lcol(2,:),'LineStyle','-','DisplayName','Zulufttemp. (MPC) [C]'); hold on;
ylabel(['$\vartheta_\mathrm{Zul}$ [$^\circ \mathrm{C}$]']); yyaxis right; title('');
plot(SimOut_PI{ix_pi}.In.T_out,'Color',lcol(1,:),'LineStyle','-','DisplayName','Außentemp. (PI) [C]'); ylabel('$\vartheta_\mathrm{Aul} [^\circ \mathrm{C}$]'); s(Nz+1).YColor = lcol(1,:); ylim([8 35]); grid on;
s(Nz+1).YColor = lcol(1,:); yyaxis left; s(Nz+1).YColor = lcol(2,:);
linkaxes(s,'x'); xlim(datetime({'2021-6-4' '2021-6-11'}))
util.formatFigure(14)
util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\MPC\LB_PiVsMpcSummer_Ts1_Hp4d.pdf','AxPosOffset',[0 -.05 0 0],'FigPosOffset',[0 0 0 -280],'SubplotYSpace',.03,'EqualX',true)

%% Plot single optimized trajectory
t = '2021-6-4';
t = '2021-12-25 12:00';
k = find(datenum(t)==datenum(getabstime(SimOut_MPC{1}.Uopt)));
lcol = lines(2); Hy = SimOut_MPC{ix_mpc}.Parameters.Ctrl.MPC.Hy; Hu = SimOut_MPC{ix_mpc}.Parameters.Ctrl.MPC.Hu;
Y = SimOut_MPC{ix_mpc}.Yopt.Data(:,:,k);
U = SimOut_MPC{ix_mpc}.Uopt.Data(:,:,k);
k = find(datenum(t)==datenum(getabstime(SimOut_MPC{ix_mpc}.In.T_ref_upper))) + SimOut_MPC{ix_mpc}.Parameters.Ctrl.MPC.T_s/SimOut_MPC{ix_mpc}.Parameters.T_s*Hy;
Ybnd = [SimOut_MPC{1}.In.T_ref_lower.Data(k)' SimOut_MPC{1}.In.T_ref_upper.Data(k)'];
k = find(datenum(t)==datenum(getabstime(SimOut_MPC{ix_mpc}.In.T_ref_upper))) + SimOut_MPC{ix_mpc}.Parameters.Ctrl.MPC.T_s/SimOut_MPC{ix_mpc}.Parameters.T_s*Hu;
Tout = SimOut_MPC{1}.In.T_out.Data(k)';
figure('Color','w','Position',[100 100 900 900]);
for nz = 1:Nz
    s(nz) = subplot(Nz+1,1,nz); stairs(Hy,Y(nz,:),'Color',lcol(2,:)); ylabel(['$\vartheta_' num2str(nz) ' [^\circ \mathrm{C}$]']); hold on; stairs(Hy,Ybnd,'LineStyle','--','Color',lcol(2,:)); hold on; 
    yyaxis right; stairs(Hu,U(nz+1,:),'Color',lcol(1,:)); ylabel(['$\dot{m}_' num2str(nz) '$ [kg/s]']); s(nz).YColor = lcol(1,:); yyaxis left; s(nz).YColor = lcol(2,:); 
    s(nz).XTickLabel = {}; xlim([0 96]);
end
s(Nz+1) = subplot(Nz+1,1,Nz+1); stairs(Hu,U(1,:),'Color',lcol(2,:)); ylabel(['$\vartheta_\mathrm{Zul}$ [$^\circ \mathrm{C}$]']); xlim([0 96]); xlabel('Pr\"{a}diktionshorizont'); hold on; stairs(Hu,[lbu(:,1) ubu(:,1)],'LineStyle','--','Color',lcol(2,:)); hold on; %ylim([14 31]); 
yyaxis right; stairs(Hu,Tout,'Color',lcol(1,:)); ylabel('$\vartheta_\mathrm{Aul} [^\circ \mathrm{C}$]'); s(Nz+1).YColor = lcol(1,:); yyaxis left; s(Nz+1).YColor = lcol(2,:);
linkaxes(s,'x'); util.formatFigure(14); 
util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\MPC\LB_MpcTrajWinter_Ts1_Hp4d.pdf','AxPosOffset',[0 -.02 0 0],'FigPosOffset',[0 0 10 -200],'SubplotYSpace',.03,'EqualX',true)


%% Plot Energy
figure('Color','w','Position',[100 100 800 300]); col = lines(3);
plot(SimOut_PI{ix_pi}.Cost_cool,'DisplayName','K\"uhlenergie (Referenz)','LineStyle','--','Color',col(1,:)); hold on;
plot(SimOut_PI{ix_pi}.Cost_heat,'DisplayName','Heizenergie (Referenz)','LineStyle','--','Color',col(2,:)); hold on;
plot(SimOut_PI{ix_pi}.Cost_el,'DisplayName','Elektroenergie (Referenz)','LineStyle','--','Color',col(3,:)); grid on;
plot(SimOut_MPC{ix_mpc}.Cost_cool,'DisplayName','K\"uhlenergie (MPC)','LineStyle','-','Color',col(1,:)); hold on;
plot(SimOut_MPC{ix_mpc}.Cost_heat,'DisplayName','Heizenergie (MPC)','LineStyle','-','Color',col(2,:)); hold on;
plot(SimOut_MPC{ix_mpc}.Cost_el,'DisplayName','Elektroenergie (MPC)','LineStyle','-','Color',col(3,:)); grid on;
title(''); ylabel('Kosten [Euro]'); legend('Location','northwest'); 
util.formatFigure(14);
set(gca,'XTickLabel',cellfun(@(x) x(1:5),get(gca,'XTickLabel'),'UniformOutput',false))
util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\MPC\LB_MpcCost_Ts1_Hp4d.pdf','AxPosOffset',[0 0 0 0],'FigPosOffset',[0 0 0 0]);

%% Plot Energy and Cost
sz = 14;
E = [SimOut_PI{ix_pi}.E_cool.Data(end)  SimOut_PI{ix_pi}.E_heat.Data(end)   SimOut_PI{ix_pi}.E_el.Data(end)];
C = [SimOut_PI{ix_pi}.Cost_cool.Data(end)    SimOut_PI{ix_pi}.Cost_heat.Data(end)     SimOut_PI{ix_pi}.Cost_el.Data(end)];
for k = 1:length(SimOut_MPC)
    E = [E; SimOut_MPC{k}.E_cool.Data(end)      SimOut_MPC{k}.E_heat.Data(end)	SimOut_MPC{k}.E_el.Data(end)];
  	C = [C; SimOut_MPC{k}.Cost_cool.Data(end)  	SimOut_MPC{k}.Cost_heat.Data(end)  	SimOut_MPC{k}.Cost_el.Data(end)];
end
xtl = {'Referenz' 'MPC (a)' 'MPC (b)' 'MPC (c)' 'MPC (d)'};
lgd = {'K\"uhlenergie' 'Heizenergie' 'Elektroenergie'};
figure('color','w','Position',[200 200 1250 330]);
subplot(1,3,1); bar(E/1000,'stacked'); ylabel('Energie [MWh]'); xticklabels(xtl); legend(lgd,'Position',[0.13,0.35,0.15,0.22]); grid on; xtickangle(45);
subplot(1,3,2); bar(C/1000,'stacked'); ylabel('Kosten [Tausend Euro]'); xticklabels(xtl); legend(lgd,'Position',[0.44,0.49,0.17,0.21]); grid on; xtickangle(45);
subplot(1,3,3);
Com = SimOut_PI{ix_pi}.Comfort_lin.Data(end); C = SimOut_PI{ix_pi}.Cost.Data(end)/1000;
scatter(Com,C,30,'bo','filled'); hold on;
text(Com+.4,C-.3,xtl{1},'Interpreter','latex','FontSize',sz);
for k = 1:length(SimOut_MPC)
    Com = SimOut_MPC{k}.Comfort_lin.Data(end,:); C = SimOut_MPC{k}.Cost.Data(end)/1000;
    scatter(Com,C,30,'ro','filled'); hold on;
    switch k
        case 1; text(Com+0.4,C-.3,xtl{k+1},'Interpreter','latex','FontSize',sz);
        case 2; text(Com-5.8,C+.2,xtl{k+1},'Interpreter','latex','FontSize',sz);
        case 3; text(Com-7,C-.3,xtl{k+1},'Interpreter','latex','FontSize',sz);
        case 4; text(Com+2,C+.2,xtl{k+1},'Interpreter','latex','FontSize',sz);
    end
end
xlabel('Komfortverletzungen [Kh]'); ylabel('Kosten [Tausend Euro]'); grid on; xlim([510 530]); ylim([10 15]);
util.formatFigure(16);
util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\MPC\LB_MpcBarCost_Ts1_Hp4d.pdf','AxPosOffset',[-0.08 0.15 0 0],'FigPosOffset',[0 0 -180 50]);

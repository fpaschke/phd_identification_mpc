% Simulate temperature room with radiator and floor heating.
% The simulation starts on P.Tstart at 00:00. The 1.1.XXXX of the corresponding
% year XXXX s assumed to be a saturday. The year thus doesnt play any role.
% A thorough description of the model can be found in the package
% SingleRoom.mo. 

clear; clc
addpath('..\..\..\..');
addpath('..\model\')
addpath('.\resource\')
addpath('..\..\LectureBuilding\sim\subfun\sig\')
addpath('..\..\LectureBuilding\sim\subfun\misc\')
addpath(genpath('.\subfun\'))

%% Create Parameters and Input Signals for Simulation
% If PeriodicOccupancy = true then Office Room scenario (seasonal Model for 
% mpc) will be simulated, if false then ConferenceRoom scenario
P = getParameters('T_s',15,'PeriodicOccupancy',false);

%% Simulate PI Control Loop for different Heatup Times
P.CtrlFun = 'ctrl_pi';
T_pre = 6:12;
for k = 1:length(T_pre)
    P.Ctrl.PI.T_pre = T_pre(k);
    SimOut_PI{k} = simulate('SingleRoom',P);
end

%% Simulate MPC Control Loop
P.CtrlFun = 'ctrl_mpc';
P.Ctrl.MPC.UseWeatherPredictionModel = false;
SimOut_MPC{1} = simulate('SingleRoom',P);
P.Ctrl.MPC.UseWeatherPredictionModel = true;
SimOut_MPC{2} = simulate('SingleRoom',P);

%% Plot 
ix_pi = 3; 
plt_detailed = false; 

% Plot Trajectories (set ZOH interpolation to timeseries objects, if quality of lines is bad)
figure('color','w','Position',[200 200 900 1200]); 
s(1) = subplot(4+2*plt_detailed,1,1); 
L1 = plot(SimOut_PI{ix_pi}.T_air,'DisplayName','PI+ZP-Regler','LineWidth',.8); hold on;
L2 = plot(SimOut_MPC{1}.T_air,'DisplayName','MPC-Regler','LineWidth',.8); grid on;
plot(SimOut_PI{ix_pi}.In.T_ref_lower,'--k','LineWidth',.8); hold on;
plot(SimOut_PI{ix_pi}.In.T_ref_upper,'--k','LineWidth',.8); hold on;
%plot(SimOut_MPC{2}.T_air,'--g','DisplayName','Luft Temp. (MPC2) [C]','LineWidth',1); grid on;
ylabel('$\vartheta_\mathrm{R} [^{\circ}$ C]'); title(''); legend([L1 L2],'Location','NorthWest'); set(gca,'XTickLabel',{}); ylim([17.5 26.5]); 
s(2) = subplot(4+2*plt_detailed,1,2);
plot(SimOut_MPC{1}.In.N_pers,'-k','DisplayName','Raumbelegung [#Pers.]','LineWidth',.8);
ylabel('Belegung [Pers.]'); grid on; title(''); set(gca,'XTickLabel',{});
s(3) = subplot(4+2*plt_detailed,1,3);
plot(100*SimOut_PI{ix_pi}.H_flo,'DisplayName','PI-Regler','LineWidth',.8); hold on;
plot(100*SimOut_MPC{1}.H_flo,'DisplayName','MPC-Regler','LineWidth',.8); 
%plot(100*SimOut_MPC{2}.H_flo,'--g','DisplayName','Heizventilpos. (MPC2) [%]','LineWidth',.8); 
grid on; ylabel('$H$ [$\%$]'); title(''); legend('Location','NorthEast'); set(gca,'XTickLabel',{}); ylim([-5 105]);
s(4) = subplot(4+2*plt_detailed,1,4); 
plot(100*(1-SimOut_PI{ix_pi}.H_win),'DisplayName','ZP-Regler','LineWidth',.8); hold on;
plot(100*(1-SimOut_MPC{1}.H_win),'DisplayName','MPC-Regler','LineWidth',.8); 
%plot(100*(1-SimOut_MPC{2}.H_win),'--g','DisplayName','Verschattungspos. (MPC2) [%]'); 
ylabel('$H_\mathrm{F}$ [$\%$]'); legend('Location','SouthWest'); grid on; title(''); ylim([-5 105]); if plt_detailed; set(gca,'XTickLabel',{}); end
if plt_detailed
    s(5) = subplot(4+2*plt_detailed,1,4);
    plot(SimOut_PI{ix_pi}.E_heat_flo+SimOut_PI{ix_pi}.E_cool_flo,'-k','DisplayName','Heizenergie (PI) [MWh]','LineWidth',.8); hold on;
    plot(SimOut_MPC{1}.E_heat_flo+SimOut_MPC{1}.E_cool_flo,'-m','DisplayName','Energie (MPC) [MWh]','LineWidth',.8); hold on;
    %plot(SimOut_MPC{2}.E_heat_flo+SimOut_MPC{2}.E_cool_flo,'--g','DisplayName','Energie (MPC2) [MWh]','LineWidth',.8); hold on;
    ylabel('Energie [kWh]'); legend; grid on; title(''); set(gca,'XTickLabel',{});
    s(6) = subplot(4+2*plt_detailed,1,5);
    plot(SimOut_PI{ix_pi}.Comfort_lin,'-k','DisplayName','Komfortverletzungen (PI) [Kh]','LineWidth',.8); hold on;
    plot(SimOut_MPC{1}.Comfort_lin,'-m','DisplayName','Komfortverletzungen (MPC) [Kh]','LineWidth',.8); hold on;
    %plot(SimOut_MPC{2}.Comfort_lin,'--g','DisplayName','Komfortverletzungen (MPC2) [Kh]','LineWidth',-0.2); hold on;
    ylabel('Komfortverletzungen [Kh]'); legend; grid on; title('');
end
linkaxes(s,'x'); %xlim([datetime(P.Tstart) datetime(datestr(datenum(P.Tstart)+7*P.N_w))])
%xlim([datetime([P.Tstart(1:4) '-3-27']) datetime([P.Tstart(1:4) '-4-9'])]); 
xlim([datetime([P.Tstart(1:4) '-2-20']) datetime([P.Tstart(1:4) '-3-8'])]); 
util.formatFigure(14);
util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\MPC\SR_PiVsMPC.pdf','XShift',-.05,'FigPosOffset',[0 0 0 -290],'AxPosOffset',[0.055 -.04 0 0],'SubplotYspace',0.04,'EqualX',true);

% Plot Scatter
figure('color','w','Position',[200 200 450 350]);
for k = 1:length(SimOut_PI)
    Com = sum(SimOut_PI{k}.Comfort_lin.Data(end)); 
	C = SimOut_PI{k}.E_heat_flo.Data(end);
    scatter(Com,C,40,'ko','filled','DisplayName','PI'); hold on;
    if k == length(SimOut_PI); C = C+10; end
    text(Com+5,C,['$T_\mathrm{vor}=$' num2str(T_pre(k)) 'h'],'Interpreter','latex','FontSize',14);
    
end
Com = sum(SimOut_MPC{1}.Comfort_lin.Data(end));     
C = SimOut_MPC{1}.E_heat_flo.Data(end);
scatter(Com,C,40,'ro','filled','DisplayName','MPC1'); hold on;
text(Com+5,C+10,'MPC (Wetter exakt)','Interpreter','latex','FontSize',14);
Com = sum(SimOut_MPC{2}.Comfort_lin.Data(end));     
C = SimOut_MPC{2}.E_heat_flo.Data(end);
scatter(Com,C,40,'bo','filled','DisplayName','MPC2'); hold on;
text(Com+5,C-0,'MPC (Wetter pr\"ad.)','Interpreter','latex','FontSize',14);
xlabel('Komfortverletzungen [Kh]'); ylabel('Energie [kWh]'); 
%xlim([300 500]); ylim([2350 2750]);
xlim([320 600]); ylim([2280 2630]);
%util.formatFigure(16); 
%util.saveTightFigure(gcf,'\\tsclient\D\Diss\Bilder\MPC\SR_PiVsMpcYear_Office.pdf','FigPosOffset',[0 0 0 10],'AxPosOffset',[0 0 0 0]);

 
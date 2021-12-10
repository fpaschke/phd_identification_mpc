clc; clear; close all;

%% PARAMS
A = [1 -.95]; 
B = [0 .1];
Nexp = 100;
rng(0); u = kron(rand(10,1),ones(10,1));
yd = filter(B,A,u,[1]);
y = yd + .1*randn(100,Nexp);
plot([yd y(:,1)])

%% ESTIMATE
for ne = 1:Nexp
    % do est with estimation of ic
    [Ai,Bi] = idModels.alg.ls.oe_iv(y(:,ne),u,1,1,1,'EstIc',true);
    B1(:,ne) = Bi{1}; A1(:,ne) = Ai{1};
    % do est without estimation of ic
    [Ai,Bi] = idModels.alg.ls.oe_iv(y(:,ne),u,1,1,1,'EstIc',false); 
    B2(:,ne) = Bi{1}; A2(:,ne) = Ai{1};
end

%% PLOT Results
figure('color','w');
s(1) = subplot(1,2,1); 
leg = bplot([A1(2:end,:)'./A(2:end) B1(2:end,:)'./B(2:end)],'VariableNames',{'a_1' 'b_1'}); grid on;
title('Mit Schätzung der AW','Interpreter','latex'); legend(leg,'Location','northwest');
yl1 = get(gca,'YLim');
s(2) = subplot(1,2,2); bplot([A2(2:end,:)'./A(2:end) B2(2:end,:)'./B(2:end)],'VariableNames',{'a_1' 'b_1'}); grid on;
title('Ohne Schätzung der AW','Interpreter','latex');
yl2 = get(gca,'YLim'); legend(leg,'Location','northwest');
set(s,'YLim',[min(yl1(1),yl2(1)) max(yl1(2),yl2(2))]);
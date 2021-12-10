clc; close all; clear;
addpath('..\..\..')

%% Params
% Param
rng(0);
N = 1e4;
Np = 20;
A = [1 -0.9];
B = [0 5*(1+A(2))];
Hp = [1 10];

% GenData and Sim
% u = kron(round(10*rand(Np,1))/10,ones(N/Np,1));
u = randn(N,1);
x = filter(B,A,u);
y = round(x);

%Plot
figure('color','w');
subplot(2,1,1); plot(y,'k'); hold on; plot(x,'--k'); hold on; 
subplot(2,1,2); plot(u,'-k');
formatFigure(11,'$t$',{'' '$u[t]$'},[],[],{{'$y[t]$' '$x[t]$'} {'$u[t]$'}});

%% Setup Model and Draw Costfunc
M_arx = idModels.NsfPolyModel(1,1,1);
M_oe = idModels.NsfPolyModel(0,1,1,0,0,1);

a = A(2)+linspace(-.02,.06,21);
b = B(2)+linspace(-.1,.1,21);
for row = 1:length(a)
    row
    M_arx.A.val = [1 a(row)];
    M_oe.E.val = [1 a(row)];
    for col = 1:length(b)
        M_arx.B.val = [0 b(col)];
        M_oe.B.val = [0 b(col)];
        for k = 1:length(Hp) 
            opt = M_arx.identifyOptions('Hp',Hp(k));
            V(row,col,k) = mean(M_arx.ts_resid(y,u,opt).^2,'omitnan');
        end
        V(row,col,length(Hp)+1) = mean(M_oe.ts_resid(y,u,M_oe.identifyOptions()).^2,'omitnan');
    end
end

%% PLOT
figure('color','w','Position',[200 200 1400 300]);
for k = 1:length(Hp)+1 
	subplot(1,length(Hp)+1,k); f = contourf(a,b,V(:,:,k)'); 
    map = colormap(flipud(gray)); colorbar(); hold on; 
    V_min = min(min(V(:,:,k)));
    [r,c] = find(V(:,:,k)==V_min);
    a_min = a(r); b_min = b(c);
    plot(A(2),B(2),'ko','MarkerSize',8); grid on;
    plot(a_min,b_min,'kx','MarkerSize',8); grid on;
    %caxis([0.1 0.3])
    %set(gca,'Colorscale','log');
end
formatFigure(17,{'$a^{[1]}$' '$a^{[1]}$' '$f^{[1]}$'},{'$b_1$' [] []},[],[arrayfun(@(k) ['ARX ($k = ' num2str(k) '$)'],Hp,'UniformOutput',0) 'OE ($k=1$)']);
%util.saveTightFigure(gcf,'D:\Diss\Bilder\Identifikation\FirstOrderQuant_Cost.pdf','AxPosOffset',[-.045 .11 0 0],'FigPosOffset',[0 0 -230 50],'SubplotXSpace',.02)

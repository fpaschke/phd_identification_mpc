clc; close all; clear;
addpath('..\..\..')

%% Params
rng(0);
sigma = 1;
N = 1e3;
Np = N;
A = [1 -0.9];
B = [0 1];
Hp = [1 5 25];

% GenData and Sim
e = sigma*randn(N,1); 
u = 1*kron(rand(Np,1),ones(N/Np,1));
y = filter(B,A,u) + filter(1,A,e);

% Setup Model and Estim CovP
for nm = 1:length(Hp)+1
    if nm <= length(Hp) % ARX
        M{nm} = idModels.NsfPolyModel(1,1,1);
        M{nm}.A.val = A;
        k = Hp(nm); 
        Lambda = sigma^2;
    else % OE
        M{nm} = idModels.NsfPolyModel(0,1,1,0,0,1);
        M{nm}.E.val = A;
        k = 1;
        Lambda = var(M{nm}.ts_resid(y,u,M{nm}.identifyOptions('Hp',1)),'omitnan');
    end
    M{nm}.B.val = B; 
    M{nm}.NoiseVariance = Lambda;
    M{nm}.Info.OptionsUsed = M{nm}.identifyOptions('Hp',k);
	M{nm}.Info.CovP = M{nm}.calcCovP(y,u,M{nm}.Info.OptionsUsed);
end

%% Show Bode Plots
M{1}.showBode('ShowFrequencyWeighting',true,'ShowH',false,'ShowPhase',false,'Models',M(2:end-1),...
    'legend',[arrayfun(@(k) ['$k = ' num2str(k) '$'],Hp,'UniformOutput',0)],'Confidence',.99)
set(gcf,'Position',[200 200 1000 350]);
C = get(gcf,'Children');
C(1).Location = 'southwest';
C(3).Location = 'southwest';
wN = pi/M{1}.Ts; 
w = linspace(floor(log10(wN/1000)),wN,1e3); % Angular Nyquist freq in 1/obj.TimeUnit
[~,~,G] = M{1}.calcConfBode(w);
s(1) = subplot(1,2,1); hold on; L = semilogx(w,20*log10(G{1}),'-k'); L.Annotation.LegendInformation.IconDisplayStyle = 'off'; hold on; 
s(1).XTickLabel = {'$10^{-3}$' '$10^{-2}$' '$10^{-1}$' '$10^{0}$'}; ylim([-20 25]); xlabel('$\omega$');
s(2) = subplot(1,2,2); s(2).XTickLabel = {'$10^{-3}$' '$10^{-2}$' '$10^{-1}$' '$10^{0}$'}; xlabel('$\omega$');
s(1).XTickLabel = {'$10^{-3}$' '$10^{-2}$' '$10^{-1}$' '$10^{0}$'};
util.formatFigure(13)
%util.saveTightFigure(gcf,'D:\Diss\Bilder\Identifikation\ARX_BodeMag.pdf','AxPosOffset',[-.07 .04 0 0],'FigPosOffset',[0 0 -120 20])

%% Draw Costfunc
clear V;
a = A(2)+linspace(-.09,.05,31);
a1 = A(2)+linspace(-.07,.07,31);
b = B(2)+linspace(-.9,.9,21);
for nm = 1:length(Hp)+1
    for row = 1:length(a)
        row
        if nm <= 3
            M{nm}.A.val = [1 a(row)];
            k = Hp(nm);
        else
            M{nm}.E.val = [1 a(row)];
            k = 1;
        end
        for col = 1:length(b)
            M{nm}.B.val = [0 b(col)];
            V(row,col,nm) = mean(M{nm}.ts_resid(y,u,M{1}.identifyOptions('Hp',k)).^2,'omitnan');
        end
    end
end

%% PLOT
figure('color','w','Position',[200 200 800 550]);
for nm = 1:length(Hp)+1 
    if nm <= 3; k = Hp(nm); else; k = 1; end
	subplot((length(Hp)+1)/2,(length(Hp)+1)/2,nm); 
	f = contourf(a,b,log10(V(:,:,nm)'));
    map = colormap(flipud(gray)); colorbar(); hold on; 
    plot(A(2),B(2),'ko'); hold on;
    M{nm}.identify(y,u,'Hp',k);
    if nm <= 3
        plot(M{nm}.A.val(2),M{nm}.B.val(2),'kx'); hold on; grid on;
    else
        plot(M{nm}.E.val(2),M{nm}.B.val(2),'kx'); hold on; grid on;
    end
end
util.formatFigure(13,{'$a^{[1]}$' '$a^{[1]}$' '$a^{[1]}$' '$f^{[1]}$'},{'$b^{[1]}$' '$b^{[1]}$' '$b^{[1]}$' '$b^{[1]}$'} ,[],[arrayfun(@(k) ['ARX ($k = ' num2str(k) '$)'],Hp,'UniformOutput',0) 'OE ($k=1$)']);
%util.saveTightFigure(gcf,'D:\Diss\Bilder\Identifikation\ARX_Cost.pdf','AxPosOffset',[-.05 -.02 0 0],'FigPosOffset',[0 0 -100 -30])

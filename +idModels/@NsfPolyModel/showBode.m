function fig = showBode(obj,varargin)
%SHOWBODE Show Bode Magnitude plot of LTI-part (G and H) of NsfPolyModel:
% y[t] = G(q^-1) g(u[t],y[t]) + H(q^-1) e[t]
%
% f = showBode(obj,varargin)
% obj [NsfPolyModel]: The Model object.
% f [Figure]:                       figure handle
% varargin:                         optional name value pairs
% omega [pos. double vector]:       Frequency points to be plotted. Def.: Automatically creates frequecy Vector
% ShowH [log.]:                     Plots Bode Magnitude plot NoiseModel H. Def. false
% ShowFrequencyWeighting [log.]:    Plots frequency weighting function H^[0,k-1]/H. Def. false
% ShowConfidence [log.]:            Plots Confidence of Magnitude and phase if true. Def. true
% Confidence [double 0...1]:        Confidence bound of bode Magnitude plots. .99 plot 99% confidence interval. Def. .99. 
%                                   .99 means that 99% of estimates will lie within that bound
% HighlightFrequencies [double]:    Highlights specific Frequencies with a vertical line in 1/obj.TimeUnit. Def.: []
% ShowPhase [0/1]:                  Plots Phase if true. Def. false.
% Labels [char]:                    Defines labeling Method. Def.: 'Full'
%                                       'Full':     Plots all labels to y axis
%                                       'Sparse':   Reduced labeling.

wN = pi/obj.Ts; % Angular Nyquist freq in 1/obj.TimeUnit
wL = floor(log10(wN/1000)); % Lower bound of frequency plot = 10^wL

p = inputParser;
addRequired(p,'obj');
addParameter(p,'omega',logspace(wL,log10(wN),1e3),@(x) isvector(x) && all(x>=0) && all(isreal(x)));
addParameter(p,'ShowFrequencyWeighting',false);
addParameter(p,'ShowH',true);
addParameter(p,'ShowConfidence',true);
addParameter(p,'Confidence',.99)
addParameter(p,'HighlightFrequencies',[],@(x) isvector(x) && all(x>0));
addParameter(p,'ShowPhase',true,@(x) x==0 || x==1);
addParameter(p,'Labels','Sparse',@(x) any(strcmpi(x,{'Full' 'Sparse'})));
addParameter(p,'Axes',[]);
addParameter(p,'Idx',[],@(x) isvector(x) && all(mod(x,1)==0));

% Undocumented yet
addParameter(p,'Models',{},@(x) iscell(x) || isa(x,'idModels.NsfPolyModel'));
addParameter(p,'YLim',[],@(x) length(x)==2 && isvector(x));
addParameter(p,'Legend',[])
parse(p,obj,varargin{:});

w = p.Results.omega;
ny = obj.OutputDimension;
nv = size(obj.B,2);
addMod = p.Results.Models; 
if ~iscell(addMod); addMod = {addMod}; end
models = [{obj} addMod];
if any(strcmpi(p.UsingDefaults,'Legend')) % Use Names if Nothing supplied
    legs = cellfun(@(m) m.Name,models,'UniformOutput',false);
else
    if isempty(p.Results.Legend)
        legs = cellfun(@(m) '',models,'UniformOutput',false); 
    else
        legs = p.Results.Legend;
    end
end

%assert(length(models) == 1 || p.Results.ShowPhase == 0,'Multiple Models cant be plotted if ShowPhase = true!');
if length(models)>1
    colors = lines(length(models));
elseif length(models)==1 && p.Results.ShowPhase
    colors = lines(2);
else
    colors = [0 0 0];
end
ShowW = p.Results.ShowFrequencyWeighting;
ShowH =  p.Results.ShowH;
ncol = nv+ShowH+ShowW;

if isempty(p.Results.Axes) 
    fig = figure('color','w');
    rows = 1:ny;
    cols = 1:ncol;
else
    subplot(p.Results.Axes);
    rows = p.Results.Idx(1);
    assert(~isempty(p.Results.Idx),'If Axes is specified the Idx need to be supplied as well!');
    if length(p.Results.Idx) == 1
        cols = nv + 1;
    elseif length(p.Results.Idx) == 2
       	cols = p.Results.Idx(2);
    else
        error('Idx needs to be a scalar or a vector!');
    end
end

for nm = 1:length(models)
    % Plot
    [confBodeMagG, confBodeMagH, G, H, confPhaseG, confPhaseH] = models{nm}.calcConfBode(w,p.Results.Confidence);
    np = 1;
    
    for row = rows
        ShowH = true;
        ShowW = false;
        for col = cols
            if isempty(p.Results.Axes) 
                s = subplot((1+(p.Results.ShowPhase)*(length(models)>1))*ny,nv+p.Results.ShowFrequencyWeighting+p.Results.ShowH,np);
            end
            % Select 
            if col <= nv && nv>0
                Amp = abs(G{row,col});
                ConfAmp = confBodeMagG{row,col};
                Phase = 180/pi*unwrap(angle(G{row,col})); 
                ConfPhase = confPhaseG{row,col};
               	if nv == 1 && ny == 1
                    lab = '$G(\omega)$'; 
                else
                    lab = ['$G_{' num2str(row) num2str(col) '}(\omega)$'];
                end
            elseif p.Results.ShowH && ShowH
                ShowH = false;
                Amp = abs(H{row,1});
                ConfAmp = confBodeMagH{row,1};
                Phase = 180/pi*unwrap(angle(H{row,1})); 
                ConfPhase = confPhaseH{row,1};
               	if ny == 1
                    lab = '$H(\omega)$'; 
                else
                    lab = ['$H_{' num2str(row) num2str(row) '}(\omega)$'];
                end
            elseif p.Results.ShowFrequencyWeighting
                ShowW = true;
            	W = models{nm}.calcFreqWeight(w);
                Amp = abs(W);
                Phase = 180/pi*unwrap(angle(W)); 
                if ny == 1
                 	lab = '$H^{[0,k-1]}(\omega)/H(\omega)$'; 
                else
                	lab =['$H^{[0,k-1]}_{' num2str(row) num2str(row) '}(\omega)/H_{' num2str(row) num2str(row) '}(\omega)$'];
                end
            end
            % Plot Amp
            semilogx(w,20*log10(Amp),'LineStyle','-','Color',colors(nm,:),'DisplayName',legs{nm}); hold on;
            % PlotConfAmp
            if p.Results.ShowConfidence && ~ShowW
                hp = semilogx(w,20*log10(Amp) + [ConfAmp -ConfAmp],'LineStyle','--','Color',colors(nm,:)); hold on;
                arrayfun(@(p) set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'),hp); % disable legend entries
            end
            %Plot Highlighted Freqencies
            for nl = 1:length(p.Results.HighlightFrequencies)
               hp = xline(p.Results.HighlightFrequencies(nl),'-y'); hold on;
               arrayfun(@(p) set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'),hp); % disabloe legend entries
            end
         	% Ylabel only on first left ax
            if col == 1 
               ylabel('Amplitude $A(\omega)$ in dB','Interpreter','latex');
            end 
            % Legend on if not empty legend and multiple models
          	if length(models)>1 && all(~cellfun(@isempty,legs))
                legend('show','Interpreter','latex');
            end
            title(lab,'Interpreter','latex'); grid on;
            xlim([w(1) w(end)]);
            s.YAxis(1).TickLabelInterpreter = 'latex';
         	if row == rows(end) && nm==1
                s.XAxis.TickLabelInterpreter = 'latex';
                xlabel(['$\omega$ in ' util.convTimeUnitsToPlotLabel(obj.TimeUnit) '$^{-1}$'],'Interpreter','latex');
            else
                s.XAxis.TickLabel = {};
            end
            
            if p.Results.ShowPhase 
                if length(models)==1
                    yyaxis right; Col = colors(2,:); 
                else
                    s = subplot((1+(length(models)>1))*ny,nv+p.Results.ShowFrequencyWeighting+p.Results.ShowH,np+nv+p.Results.ShowFrequencyWeighting+p.Results.ShowH);
                    Col = colors(nm,:); 
                end
                semilogx(w,Phase,'LineStyle','-','Color',Col,'DisplayName',legs{nm}); hold on; grid on;
                if p.Results.ShowConfidence && ~ShowW
                    hp = semilogx(w,Phase + 180/pi*[ConfPhase -ConfPhase],'LineStyle','--','Color',Col); hold on;
                    arrayfun(@(p) set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'),hp); % disable legend entries
                end
                if length(models)==1
                    set(gca,'YColor',colors(nm+1,:));
                    s.YAxis(1).Color  = colors(1,:);
                    s.YAxis(2).TickLabelInterpreter = 'latex';
                    if col == cols(end) 
                        ylabel('Phase $\varphi(\omega)$ [$^{\circ}$]','Interpreter','latex'); 
                    end
                elseif length(models)>1 
                    if col == 1
                        ylabel('Phase $\varphi(\omega)$ [$^{\circ}$]','Interpreter','latex');
                    end
                    if all(~cellfun(@isempty,legs))
                        legend('show','Interpreter','latex');
                    end
                    xlim([w(1) w(end)]);
                	if row == rows(end)
                        s.XAxis.TickLabelInterpreter = 'latex';
                        xlabel(['$\omega$ in ' util.convTimeUnitsToPlotLabel(obj.TimeUnit) '$^{-1}$'],'Interpreter','latex');
                    else
                        s.XAxis.TickLabel = {};
                    end
                end
                s.YAxis(1).TickLabelInterpreter = 'latex';
            end   
            if length(s.YAxis)>1; s.YAxis(2).TickLabelInterpreter = 'latex'; end 
            np = np + 1;
        end
    end  
end

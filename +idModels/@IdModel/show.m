function h = show(obj,type,y,varargin)
%Plots residual statistics specified by type. Type can be one of the following:
% 'Autocovariance': Plots unbiased estimate of autocovariance function of prediction errors.
% 'Boxplot': Plots boxplot of prediction errors. Mean will be plotted with +.
% 'Histogram': Plots histogramm of prediction errors.
% 'Scatter': Plots Predicted Model Output vs Measured values. 
% 'Prediction': Plots Timeseries of predicted vs measured values. Each Dataset will be plotted in seperate TAB.
%In case of dynamic model the Hp-Step-ahead prediction errors will be shown. Creates a subplot for each output.
%
% show(obj,type,y,varargin)
% obj [BModel]: Any identified Model with available simulate routine
% y [(cell of) N x ny double or signals.TSignalCollection]: Measured OutputData or TSignalCollection Object containing u and y Data.
% u [(cell of) N x nu double]: Measured InputData
% varargin: Name value pairs -> see Below
%
% Opt. Name Value Pairs that are valid for all Plot types
% Axes [Axes (handle)]: Axes to which the plots will be drawn (not available for type = 'Prediction').
% Hp [pos int]: Prediction Horizon in #Samples. Only available for dynamic models. (Default 1)
% ObserverInitSamples [pos int]: #Samples to compute initial value of observer. Only available for dynamic models. (Default 'auto')
% FontSize [pos. double]: Determines Fontsize used for axes and legends (Def. 12)
% IxOut [pos. int vector]: Defines indices of outputs to be plotted  
%
% Opt. Name Value Pairs for type = 'Boxplot'
% Color [char]: Determines Color of lines (Def. 'k')
% ShowExpStd [char]: Plots expected Standard Dev. of the Model Errors with specified LineSpec (e. g. '--r'). [] does not plot Std. Dev.
%                    (Only available for dyn. models without output feedback) (Def. []).
% Whisker [double in range of 0 to 100]: Determines position of whiskers in %. 2.5 means upper and lower 2.5% quantile, i. e.
%                                        95% of the prediction errors will be within the whiskers. (Def. 5) 
% ShowOutliers [logical]: Plots markers for outliers if true. (Def. false)
% DeltaK [pos int]: Plots boxes only that are DeltaK steps away. 2 will plot every second box only. (Def.: 1) 
%
% Opt. Name Value Pairs for type = 'Histogramm'
% Bins [pos. int]: Number of Bins of Boxplot. Def. 50
% ShowNormProb [logical]: Adds probability values corresponding to normal distribution if true (Def. false).
%
% Opt. Name Value Pairs for type = 'Autocovariance'
% Maxlag [pos. int]:    Maximum lag of Autocovariance function. (Def. 25)
% Color [char]:         Determines Color of lines (Def. 'k')
% Significance [double]:Statistical significance of critical values of autocovariance r(tau) (0< alpha < 1) [Def. .01]. 
%                       Under the assumption of residual whiteness  (Null hypothesis), r(tau) will be outside of r_crit
%                       with probability alpha (higher values of alpha lead to lower r_crit)
% LineStyle [char]:     Any Valid Linemarker e.g. '-' , '--', '-.'. [Def.: '-']                              
%
% Opt. Name Value Pairs for type = 'Scatter'
% Alpha [pos. int]: Transperency value for Markers 1=no transperancy. (Def. 0.05)
%
% Opt. Name Value Pairs for type = 'Prediction'
% ShowInputs [cell of pos int vectors]: Determines appearance of plot: {[1 3] [2]} plots 2 extra subplots 1st containing
%                                       Input number 1 and 3 and second with input 2. [Def. {{1} {2} ... {obj.InputDimension}}]
% LeftAx [cell of pos logical vectors]: Determines appearance of plot: {[1 0] [1]} plots 2 extra subplots 1st containing
%                                       Input number 1 on left ax and input 3 on right ax and second with input 2 on left ax.
%                                       [] plots everything on left axes. (Def. [])
% LineSpec [cell of cellstr]:   Determines plotmarkers of each input: {{'-k' '--r'} {-m}} plots line 1 and 2 in subplot 1 with 
%                               black and red dashed line resp. [] uses matlab defaults. (Def. [])
% HighlightOutliersByNsigma [pos. double]:  Highlights Outliers with 'ro'. If HighlightOutliersByNsigma=3 then an outlier is
%                                           defined by 3*(empirical) Standarddeviation. [] means no highlighting. (Def. [])
% TimeUnit ['char']: Set specific TimeUnit of x-axis. (Def. obj.TimeUnit)
%
% EXAMPLES:
% 1.) Softplus Model
% M = idModels.SoftplusModel('Name','SoftPlus','Parameters',[2 1 10],'Free',[1 1 1]);
% u = (-50:.1:50)';
% rng(0); y = 3*(u-20).*(u>=20) + 100 + 5*randn(length(u),1);
% M.identify(y,u,'EstimateOutputOffset',true);
% M.showModel('Data',{y u})
% M.show('Histogram',y,u,'Bins',50,'ShowNormProb',1)
% M.show('Prediction',y,u,'HighlightOutliersByNsigma',2)
% M.show('Boxplot',y,u)
% M.show('AutoCovariance',y,u)
% M.show('Scatter',y,u)
%
% 2.) Dynamic MISO-ARMAX with 3 Inputs
% u = kron(rand(30,3),ones(10,1)); e = .1*randn(length(u),1);  
% y = filter([0 .1],[1 -.9],u(:,1)) + filter([0 -.1],[1 -.9],u(:,2)) + filter([0 .05],[1 -.9],u(:,3)) + filter([1 .8],[1 -.9],e); % ARMAX of order 1
% u = {u(1:200,:) u(201:end,:)}; y = {y(1:200,:) y(201:end,:)}; % Split for Demo
% M = idModels.NsfPolyModel(1,[1 1 1],[1 1 1],1);
% M.identify(y,u);
% M.show('Boxplot',y,u,'Hp',10,'ShowExpStd','--r')
% M.show('Prediction',y,u,'Hp',2,'ShowInputs',{[1 2] [3]},'LeftAx',{[0 1] [0]},'LineSpec',{{'-r' '-b'} {'--k'}},'HighlightOutliersByNsigma',2);

%% PREPARE DATA AND INIT
p = inputParser;
addRequired(p,'Type',@(x) any(strcmpi(x,{'Autocovariance' 'Boxplot' 'Histogram' 'Scatter' 'Prediction'})));
addRequired(p,'y',@(x) isa(x,'signals.TSignalCollection') || isnumeric(x) || all(cellfun(@isnumeric,x)));
addOptional(p,'u',[],@(x) isnumeric(x) || all(cellfun(@isnumeric,x)));

% All
addParameter(p,'Hp',1,@(x) x>0 && mod(x,1)==0);
addParameter(p,'ObserverInitSamples','auto',@(x) (x>0 && mod(x,1)==0) || strcmpi(x,'auto') || isinf(x));
addParameter(p,'ObserverInitUseGls',true,@(x) x==1 || x==0);
addParameter(p,'Axes',[]);
addParameter(p,'FontSize',12,@(x) isscalar(x) && x > 1);
addParameter(p,'Color','k',@(x) ischar(x) || (isvector(x) && length(x)==3)); %boxplot and autocovariance

% Boxplot
addParameter(p,'ShowExpStd',[],@(x) islinespec(x));
addParameter(p,'Whisker',5,@(x) isscalar(x) && x>0 && x<100);
addParameter(p,'ShowOutliers',false,@(x) x==0 && x==1);
addParameter(p,'DeltaK',1,@(x) mod(x,1)==0 && x>0);
addParameter(p,'IxOut',1:obj.OutputDimension,@(x) all(arrayfun(@(ni) any(ni==1:obj.OutputDimension),x))); % Boxplot, Histogramm, AutoCov

% Histogramm 
addParameter(p,'Bins',50);
addParameter(p,'ShowNormProb',false,@(x) x==1 || x==0);

% Autocovariance 
addParameter(p,'Maxlag',25,@(x) x>0 && mod(x,1)==0);
addParameter(p,'Normalization','unbiased',@(x) strcmpi(x,{'none' 'unbiased' 'biased' 'coeff'}));
addParameter(p,'Significance',.01,@(x) x>0 && x<1);
addParameter(p,'LineStyle','-')

% Scatter 
addParameter(p,'Alpha',.05,@(x) isscalar(x) && x>=0 && x<=1);

% Prediction
addParameter(p,'ShowInputs',num2cell(1:obj.InputDimension))
addParameter(p,'LeftAx',[])
addParameter(p,'LineSpec',[])
addParameter(p,'HighlightOutliersByNsigma',[],@(x) isempty(x) || x > 0  )
if isa(obj,'idModels.DynamicModel'); addParameter(p,'TimeUnit',obj.TimeUnit); end
addParameter(p,'ShowDrift',false,@(x) x==0 && x==1);

% Parse Inputs
parse(p,type,y,varargin{:});

% Checks
if ~isempty(p.Results.Axes)
    assert(length(p.Results.IxOut) == length(p.Results.Axes),'Number of supplied Axes Objects needs to match length of IxOut (if supplied) or OutputDimension of the model.'); 
end

%% Get Data 
Hp = p.Results.Hp;
[y,u] = obj.getRawData(y,p.Results.u); % Data will be checked in calcResiduals
if isa(obj,'idModels.StaticModel') 
    [e,yp,misc] = obj.calcResiduals(y,u);
    if ~any(strcmpi(p.UsingDefaults,'Hp')) || ~any(strcmpi(p.UsingDefaults,'ObserverInitSamples'))
        warning('The Parameters "Hp" and "ObserverInitSamples" will be ignored because they are only available for dynamic models!');
        Hp = 1;
    end
else
    [e,yp,misc] = obj.calcResiduals(y,u,Hp,p.Results.ObserverInitSamples,~p.Results.ObserverInitUseGls);
end
        
%% Now make plots
if isempty(p.Results.Axes); h = figure('color','w'); end
if  strcmpi(type,'prediction')
    assert(any(strcmpi(p.UsingDefaults,'IxOut')),'Unsupported yet!');
    % Extract options
    ins = p.Results.ShowInputs; lax = p.Results.LeftAx;
    if ~iscell(ins); ins = mat2cell(ins(:),ones(1,length(ins)),1); end
    if isempty(lax); lax = cellfun(@(in) true(1,length(in)) ,ins,'UniformOutput',false); end
    if isa(obj,'idModels.DynamicModel')
        tu = p.Results.TimeUnit;
        abb = [idModels.util.getSupportedTimeUnits;
               {'ns' '$\mue$s' 'ms' 's' 'min' 'h' 'd' 'w' 'm' 'y'}];
        tu = abb{2,strcmpi(abb(1,:),tu)};
    else
        tu = '';
    end
    
    % Do plots
    if length(y) > 1; tg = uitabgroup; end
    for ns = 1:length(y)
        if length(y)>1; thistab = uitab(tg,'Title',['Set' num2str(ns)],'BackgroundColor','w'); end
        Nsp = obj.OutputDimension+length(ins);
        if isa(obj,'idModels.DynamicModel')
            t = (0:obj.Ts:obj.Ts*(length(y{ns}) -1))*util.getTimeUnitConvFactor(obj.TimeUnit)/util.getTimeUnitConvFactor(p.Results.TimeUnit);
        else
            t = 0:length(y{ns}) -1;
        end
        for nsub = 1:Nsp
            if length(y)>1; s(nsub) = subplot(Nsp,1,nsub,'Parent',thistab); else; s(nsub) = subplot(Nsp,1,nsub); end
            if nsub <= obj.OutputDimension % Plot Outputs and Predicted output
                stairs(t,y{ns}(:,nsub),'DisplayName',getOutLab(obj,nsub)); hold on;
                l = plot(t,yp{ns}(:,nsub,Hp),'DisplayName',['Pr\"adiktion ' getOutLab(obj,nsub,obj.Name)]);
                if ~isequal(p.Results.Color,'k') && ~isequal(p.Results.Color,zeros(1,3))  
                    l.Color = p.Results.Color;
                end
                if ~isempty(p.Results.HighlightOutliersByNsigma)
                    ix = abs(e{ns}(:,nsub,Hp)) > p.Results.HighlightOutliersByNsigma*misc.Std(Hp,nsub);
                    plot(t(ix),yp{ns}(ix,nsub,Hp),'ro','DisplayName','Ausrei{\ss}er');
                end
            else % Plot Inputs
                nsub_in = nsub-obj.OutputDimension;
                for nl = 1:length(ins{nsub_in})
                    if any(lax{nsub_in}~=1)
                        if lax{nsub_in}(nl) == 0; yyaxis right; else; yyaxis left; end
                    end

                    if ~isempty(p.Results.LineSpec)
                        l = stairs(t,u{ns}(:,ins{nsub_in}(nl)),p.Results.LineSpec{nsub_in}{nl},'DisplayName',getInLab(obj,ins{nsub_in}(nl))); hold on; 
                    else
                        l = stairs(t,u{ns}(:,ins{nsub_in}(nl)),'DisplayName',getInLab(obj,ins{nsub_in}(nl))); hold on;
                    end
                    if any(lax{nsub_in}==0) % if some line on right
                        if sum(lax{nsub_in}==1) == 1 && lax{nsub_in}(nl) == 1 % Only one plot on left and current is left
                            s(nsub).YAxis(1).Color  = l.Color;
                        elseif sum(lax{nsub_in}==0) == 1 && lax{nsub_in}(nl) == 0 % Only one plot on right and current is right
                            s(nsub).YAxis(2).Color  = l.Color;
                        end
                    else %black
                        s(nsub).YAxis(1).Color = 'k'; if length(s(nsub).YAxis) == 2; s(nsub).YAxis(2).Color = 'k'; end 
                    end
                end
            end
            lgd = legend('show'); set(lgd,'Interpreter','latex','Location','NorthEast','FontSize',p.Results.FontSize);
            if nsub == Nsp
                xl = 'Zeit $t$ in '; 
                if ~isempty(tu); xl = [xl tu]; end
            else
                xl = []; 
            end
            util.formatFigure(p.Results.FontSize,xl,[],[],[],[],'Axes',s(nsub));
            if nsub ~= Nsp s(nsub).XTickLabel = ''; end
        end
        linkaxes(s,'x'); xlim([0 t(end)]);
    end
else
    if iscell(e); e = cell2mat(e); end
    if iscell(yp); yp = cell2mat(yp); end
    ym = cell2mat(y);
    ncol = ceil(sqrt(length(p.Results.IxOut)));
    nrow = ceil(length(p.Results.IxOut)/ncol);

    for no =  p.Results.IxOut
        % Create Subplot each output
        if isempty(p.Results.Axes); s = subplot(nrow,ncol,no); else; s = subplot(p.Results.Axes); end
        if obj.OutputDimension > 1; lab = ['\varepsilon_{' num2str(no) '}']; else; lab = '\varepsilon'; end
        if isa(obj,'idModels.DynamicModel'); lab = [lab '[t|t-' num2str(Hp) ']']; end
        switch lower(p.Results.Type)
            case 'boxplot'
                if p.Results.ShowOutliers; so = 'outliers'; else; so = 'nooutliers'; end
                ix = p.Results.DeltaK:p.Results.DeltaK:Hp;
                util.bplot(squeeze(e(:,no,ix)),so,'whisker',p.Results.Whisker,'colors',p.Results.Color,'LineStyle',p.Results.LineStyle); 
                hold on; xlim([0 length(ix)+1]);
                set(gca,'XTick',linspace(0,length(ix),4),'XTickLabel',arrayfun(@(x) num2str(p.Results.DeltaK*x),linspace(0,length(ix),4),'UniformOutput',0));
%              	if isa(obj,'idModels.DynamicModel')
%                     ax = gca; ax.XTickLabel = arrayfun(@(xi) num2str(xi),ax.XTick*obj.Ts*util.getTimeUnitConvFactor(obj.TimeUnit)/util.getTimeUnitConvFactor(p.Results.TimeUnit),'UniformOutput',0);
%                 end
                if isa(obj,'idModels.DynamicModel'); xl = 'Pr\"adiktionshorizont $k$'; else; xl = ''; end
                util.formatFigure(p.Results.FontSize,xl,['Pr\"adiktionsfehler $' strrep(lab,['t-' num2str(Hp)],'t-k') '$'],[],[],[],'Axes',s);
                % Plot std 
                if isa(obj,'idModels.DynamicModel') 
                    plot(misc.Std(ix,no),'DisplayName','Empirische Standardabweichung','LineStyle','-','Color',p.Results.Color); hold on; 
                    if ~isempty(p.Results.ShowExpStd) 
                        if isa(obj,'idModels.NsfPolyModel') 
                            assert(isdiag(double(obj.NoiseVariance>1e-10)) && ~any(obj.HasOutputFeedback),'"NoiseVarance" isnt a diag Matrix! -> Unsupported yet!');
                            F = obj.calcKStepPredictor(Hp); Ff = F{no,no};
                            Std_exp = sqrt(obj.NoiseVariance(no,no).*cumsum((Ff.^2),2))';
                            plot(Std_exp(ix),p.Results.ShowExpStd,'DisplayName','Erwartete Standardabweichung','LineStyle',p.Results.LineStyle); hold on;      
                        else
                            error('Not Implemented yet!');
                        end
                    end
                    lgd = legend('show'); set(lgd,'Interpreter','latex','Location','NorthWest','FontSize',p.Results.FontSize);
                end
            case 'histogram'
              	mue = misc.Mean(Hp,:); stddev = misc.Std(Hp,:);
                Bins = cell2mat(arrayfun(@ (no) linspace(mue(no) - 3*stddev(no),mue(no) + 3*stddev(no),p.Results.Bins),(1:obj.OutputDimension)','UniformOutput',false));
                if min(e(:,no,Hp)) < Bins(1); Bins = [min(e(:,no,Hp)) Bins]; end
                if max(e(:,no,Hp)) > Bins(end); Bins = [Bins max(e(:,no,Hp))]; end
                histogram(e(:,no,Hp),Bins(no,:),'FaceColor',ones(1,3)*0.75,'Normalization','probability'); xlim([Bins(no,2)-stddev(no)/2  Bins(no,end-1)+stddev(no)/2]); hold on; 
                util.formatFigure(p.Results.FontSize,['Pr\"adiktionsfehler $' lab '$'],'Relative H\"aufigkeit',[],[],[],'Axes',s);
                if p.Results.ShowNormProb
                    P = diff(normcdf(Bins(no,2:end-1),mue(no),stddev(no)));
                    lin(1) = plot((Bins(no,3:end-1)+Bins(no,2:end-2))/2,P,'-k','DisplayName','Normalverteilung'); hold on;
                    if isfield(misc,'Std_exp')
                        P_exp = diff(normcdf(Bins(no,2:end-1),0,misc.Std_exp(Hp,no)));
                        lin(2) = plot((Bins(no,3:end-1)+Bins(no,2:end-2))/2,P_exp,'--k','DisplayName','Erwartete Normalverteilung'); hold on;
                    end
                    legend(lin,'Interpreter','latex','Location','NorthWest','FontSize',p.Results.FontSize);
                end
            case 'autocovariance'
                ee = e(~isnan(e(:,no,Hp)),no,Hp); 
                [r_ee,lag] = xcorr(ee,p.Results.Maxlag,'biased');
                plot(lag(p.Results.Maxlag+1:end),r_ee(p.Results.Maxlag+1:end),'LineStyle',p.Results.LineStyle,'Color',p.Results.Color); xlim([0 lag(end)]); hold on;
%                 if isa(obj,'idModels.DynamicModel')
%                     ax = gca; ax.XTickLabel = arrayfun(@(xi) num2str(xi),ax.XTick*obj.Ts*util.getTimeUnitConvFactor(obj.TimeUnit)/util.getTimeUnitConvFactor(p.Results.TimeUnit),'UniformOutput',0);
%                 end
                util.formatFigure(p.Results.FontSize,'Lag $\tau$',['Autokovarianz $\hat{r}_{\varepsilon}[\tau]$'],[],[],[],'Axes',s);
                if ~isempty(p.Results.Significance) && Hp == 1 
                    r_crit = 1./sqrt((length(ee)-(1:p.Results.Maxlag)))*util.norminv(1 - p.Results.Significance/2,0,obj.NoiseVariance(no,no));  
                    lin = plot(1:p.Results.Maxlag,r_crit,'LineStyle','-.','Color',p.Results.Color,'DisplayName',['Kritischer Wert ($\alpha=$' num2str(p.Results.Significance) ')']); hold on;
                 	plot(1:p.Results.Maxlag,-r_crit,'LineStyle','--','Color',p.Results.Color); 
                    legend(lin,'Interpreter','latex','Location','NorthWest','FontSize',p.Results.FontSize);
                end
            case 'scatter'
                scatter(ym,yp(:,no,Hp),12,'ko','filled','MarkerFaceAlpha',p.Results.Alpha,'MarkerEdgeAlpha',p.Results.Alpha,'MarkerEdgeColor',p.Results.Color,'MarkerFaceColor',p.Results.Color); hold on;
                if obj.OutputDimension == 1; yl = '$y[t]$'; else; yl = ['$y_{' num2str(no) '}[t]$']; end 
                util.formatFigure(p.Results.FontSize,['Messwert ' yl],['Pr\"adiktion $' strrep(lab,'\varepsilon','\hat{y}') '$'],[],[],[],'Axes',s);
                xl = get(s,'XLim'); yl = get(s,'YLim'); 
                plot([max(xl(1),yl(1)) min(xl(2),yl(2))],[max(xl(1),yl(1)) min(xl(2),yl(2))],'LineStyle','-','Color',p.Results.Color); hold on;
                plot([max(xl(1),yl(1)) min(xl(2),yl(2))],[max(xl(1),yl(1)) min(xl(2),yl(2))] + 3*misc.Std(Hp,no),'LineStyle','--','Color',p.Results.Color); hold on;
                lin = plot([max(xl(1),yl(1)) min(xl(2),yl(2))],[max(xl(1),yl(1)) min(xl(2),yl(2))] - 3*misc.Std(Hp,no),'LineStyle','--','Color',p.Results.Color,'DisplayName',['$\pm 3 \sigma_{' lab '}$']); ylim(yl);
                legend(lin,'Interpreter','latex','Location','NorthWest','FontSize',p.Results.FontSize);
        end
    end
end
end

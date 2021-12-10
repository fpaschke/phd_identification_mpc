function [h, sct] = showModel(obj,varargin)
%SHOW Plots the Response of the model to a separte subplot window. Function
%destinguishes two modes -> see below.
%
% h = showModel(obj,linespec,varargin)
% Plots the response of the Model wrt to each inputsignal. The inputsignal to channel i will be 
% computed by ui = linspace(obj.InputMin(i),obj.InputMax(i),Npts) where Np is the number of points
% to be set by Name Value Pair options. The other inputs will all be set to 0 for simulation except
% if Name Value Parameter U0 is supplied (The values of U0 will be used instead).
% obj [StaticModel]: Identified (no NaN Parameter values) Model to be plotted
% linespec [char]: any valid matlab linespecifier eg --k -r...
% varargin: Optional Name Value Pairs -> see below. Unmatched Parameters will be passed to axesOptions.
%
% Available Name value pairs:
% 'Npts' [pos. int scalar]: Defines the number of points to be plottes. Only available in Mode 2. Default: 50
% 'XLabel' [string]: The x axis labels. Default: obj.InputName
% 'YLabel' [string]: The y axis labels. Default: obj.OutputName
% 'Axes' [Axes handle]: Defines the axes where all the plots will be plotted in. (Def. creates nex axes)
% 'MeshPlot' [logical]: Plots 3d Meshplot if the Model has 2 Inputs. (Def. true)
% 'InputIdx' [pos. int Vector]: Defines the input channel indexes which will be plotted (Def.: 1:obj.InputDimension) (only if meshplot = 0)
% 'OutputIdx' [pos. int Vector]: Defines the output channel indexes which will be plotted (Def.: 1:obj.OutputDimension)
% 'U0' [1 x obj.InputDimension double]:  
% 'Data' [signals.TSignalCollection]:   If spüecified the rawdata values will be extracted from the Collection and will be
%                                       plotted as scatterplot (Def. []). Only Available for models with one input.

p = inputParser; p.KeepUnmatched = true;
%addOptional(p,'u',[],@(x) size(x,2)==obj.InputDimension);
addOptional(p,'LineSpec','-k',@util.islinespec)
addParameter(p,'InputIdx',1:obj.InputDimension,@(x) all(x<=obj.InputDimension & x>0 & mod(x,1)==0));
addParameter(p,'OutputIdx',1:obj.OutputDimension,@(x) all(x<=obj.OutputDimension & x>0 & mod(x,1)==0));
addParameter(p,'Npts',30,@(x) x>1 & mod(x,1)==0);
addParameter(p,'U0',zeros(1,obj.InputDimension),@(x) isvector(x) && all(isnumeric(x)));
addParameter(p,'XLabel',[],@(x) ischar(x) || iscell(x));
addParameter(p,'YLabel',[],@(x) ischar(x) || iscell(x));
addParameter(p,'ZLabel',[],@(x) ischar(x));
addParameter(p,'Axes',[]);
addParameter(p,'Data',[],@(x) isa(x,'signals.TSignalCollection') || iscell(x))
if obj.InputDimension == 2 mp = true; else mp = false; end
addParameter(p,'MeshPlot',mp,@(x) x==1 || x==0)
parse(p,varargin{:});
assert(obj.InputDimension == 1 || isempty(p.Results.Data),'Data can be used only for Models with one input!');
if isempty(p.Results.Axes)
    h = figure('color','w');     
end

sct = [];
np = 1; 
args = [util.struct2namevaluepairs(p.Unmatched) 'Interpreter','latex','Legend','off'];
for no = p.Results.OutputIdx
    flag = true;
    if p.Results.MeshPlot
        assert(obj.InputDimension == 2,'"MeshPlot" is only available for 2 Input systems!');
        if isempty(p.Results.Axes)
            subplot(length(p.Results.OutputIdx),1,np)
        else
            subplot(p.Results.Axes) 
        end
        x = linspace(obj.InputMin(1),obj.InputMax(1),p.Results.Npts)';
        y = linspace(obj.InputMin(2),obj.InputMax(2),p.Results.Npts)';
        for k = 1:length(x)
            z(:,k) = obj.simulate([x(k)*ones(length(x),1) y]);
        end
        hold on; surf(x,y,z,'DisplayName',[util.convToPlotLabel(obj.Name,'latex') ' ($\sigma_e$=' num2str(sqrt(obj.NoiseVariance)) ')']); grid on;


        if  isempty(get(get(gca,'XLabel'),'String'))
            if isempty(p.Results.XLabel)
                xl = util.convToPlotLabel(obj.InputName{1},'latex');
            else
                xl = p.Results.XLabel;
            end
            args = [args 'XLabel' xl];
        end
        if  isempty(get(get(gca,'YLabel'),'String'))
            if isempty(p.Results.YLabel)
                xl = util.convToPlotLabel(obj.InputName{2},'latex');
            else
                xl = p.Results.YLabel;
            end
            args = [args 'YLabel' xl];
        end
        if  isempty(get(get(gca,'ZLabel'),'String'))
            if isempty(p.Results.ZLabel)
                xl = util.convToPlotLabel(obj.OutputName{no},'latex');
            else
                xl = p.Results.ZLabel;
            end
            args = [args 'ZLabel' xl];
        end
        util.setAxesOptions(util.axesOptions(args{:}));
    else
        for ni = p.Results.InputIdx
            u = linspace(obj.InputMin(ni),obj.InputMax(ni),p.Results.Npts)';
            U = [p.Results.U0(1:ni-1).*ones(p.Results.Npts,ni-1) u p.Results.U0(ni+1:end).*ones(p.Results.Npts,obj.InputDimension-ni)];
            y = obj.simulate(U);
            if isempty(p.Results.Axes)
                subplot(length(p.Results.OutputIdx),length(p.Results.InputIdx),np)
            else
                figure(p.Results.Axes.Parent)
                subplot(p.Results.Axes) 
            end
            hold on;
         	if ~isempty(p.Results.Data) && obj.InputDimension == 1
                if iscell(p.Results.Data)
                    gca; 
                    yd = p.Results.Data{1}; if iscell(yd) yd = cell2mat(yd); end
                    x = p.Results.Data{2}; if iscell(x) x = cell2mat(x); end
                    sct = scatter(x,yd,12,'filled','k');
                else
                    tsc = p.Results.Data;
                    sct = tsc(obj.InputName{1}).scatter(tsc(obj.OutputName{no}),'Axes',gca);
                end
            end   
            hold on;
            if flag
                flag = false;
                plot(u,y(:,no),p.Results.LineSpec,'DisplayName',[util.convToPlotLabel(obj.Name,'latex') ' ($\sigma_e$=' num2str(sqrt(obj.NoiseVariance(no,no))) ')']);
            else
                hl = plot(u,y(:,no),p.Results.LineSpec);
                set(get(get(hl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % disable legend entries
            end
            grid on;

            if  isempty(get(get(gca,'XLabel'),'String'))
                if isempty(p.Results.XLabel)
                    xl = ['From: ' util.convToPlotLabel(obj.InputName{ni},'latex')];
                else
                    if ischar(p.Results.XLabel) xl = p.Results.XLabel; else xl = p.Results.XLabel{ni}; end;
                end
                args = [args 'XLabel' xl];
            end
            if isempty(get(get(gca,'YLabel'),'String'))
                if isempty(p.Results.YLabel)
                    yl = ['To: ' util.convToPlotLabel(obj.OutputName{no},'latex')];
                else
                    if ischar(p.Results.YLabel) yl = p.Results.YLabel; else yl = p.Results.YLabel{no}; end;
                end
                args = [args 'YLabel' yl]; 
            end
            util.setAxesOptions(util.axesOptions(args{:}))
            np = np + 1;
        end
    end
end

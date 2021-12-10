function formatFigure(sz,xlab,ylab,zlab,tit,leg,varargin)
%Format all axes of current figure.
%
% formatFigure(sz,xlab,ylab,zlab,tit,leg,varargin)
% sz [pos int]: Size of labels
% xlab [char or cell of char]: x axes labels to be set.
% ylab [char or cell of char]: y axes labels to be set.
% zlab [char or cell of char]: z axes labels to be set.
% tit [char or cell of char]: titles of each subplot.
% leg [char or cell of char]: legends of each subplot.
% varargin: Optional Name Value Pairs
%   'LegLoc' [char]: Any valid argument for 'Location' property of Legend (e. g. 'NorthWest', 'SouthEast', etc.)
%   'Axes' [Axes object]: If supplied then only this specific Axes object will be formated.

p = inputParser();
addParameter(p,'Axes',[])
addParameter(p,'LegLoc','NorthEast')
parse(p,varargin{:});

if isempty(p.Results.Axes)
    f = gcf;
    a = findall(f,'type','axes');
    a = a(end:-1:1);
else
    a = p.Results.Axes;
end
nAx = length(a);
if nargin > 1 && ~iscell(xlab); xlab = arrayfun(@(ni) xlab,1:nAx,'UniformOutput',false); end
if nargin > 2 && ~iscell(ylab); ylab = arrayfun(@(ni) ylab,1:nAx,'UniformOutput',false); end
if nargin > 3 && ~iscell(zlab); zlab = arrayfun(@(ni) zlab,1:nAx,'UniformOutput',false); end
if nargin > 4 && ~iscell(tit); tit = arrayfun(@(ni) tit,1:nAx,'UniformOutput',false); end
if nargin > 5 && (~iscell(leg) || nAx ~= length(leg)); leg = arrayfun(@(ni) leg,1:nAx,'UniformOutput',false); end 

for na = 1:nAx
    a(na).FontSize = sz;
    grid(a(na),'on');
    cb = a(na).Colorbar;
    if ~isempty(cb)
        cb.FontSize = sz;
        cb.TickLabelInterpreter = 'latex';
    end
   
  	xlabel(a(na),a(na).XLabel.String,'Interpreter','latex','FontSize',sz);
    a(na).XAxis.TickLabelInterpreter = 'latex';
    if nargin > 1; xlabel(a(na),xlab{na}); end
    
	a(na).YAxis(1).Label.Interpreter = 'latex';
    a(na).YAxis(1).TickLabelInterpreter = 'latex';
    if length(a(na).YAxis) == 2 
        a(na).YAxis(2).TickLabelInterpreter = 'latex'; 
        a(na).YAxis(2).Label.Interpreter = 'latex';
    end
    if nargin > 2; ylabel(a(na),ylab{na}); end
    
 	zlabel(a(na),a(na).ZLabel.String,'Interpreter','latex','FontSize',sz);
	a(na).ZAxis.TickLabelInterpreter = 'latex';
    if nargin > 3; zlabel(a(na),zlab{na}); end
    
    if nargin > 4
     	title(a(na),tit{na},'Interpreter','latex','FontSize',sz);
    end
    
   	if nargin > 5 && ~isempty(leg{na})
     	legend(a(na),leg{na});
    end
    if ~isempty(a(na).Legend)
        if any(strcmpi(p.UsingDefaults,'LegLoc'))
            legend(a(na),a(na).Legend.String,'Interpreter','latex','FontSize',sz);
        else
            legend(a(na),a(na).Legend.String,'Interpreter','latex','FontSize',sz,'Location',p.Results.LegLoc);
        end
    end
end

return;
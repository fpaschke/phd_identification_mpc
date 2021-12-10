function axopt = axesOptions(varargin)
%Returns available Options for 2D-Axes from name value pairs. The returned structure
%axopt contains available name value pairs and a field unmatched containing the cell of
%unmatched name value pairs.
%
% axopt = axesOptions(varargin)
%
% Avail. Name Value Pairs:
% 'Position': determines the position of the axes object. 2 possible Methods are available to specify the postion:
%       Method 1:   double vector with 4 elements indicating [x y w h] where x and y
%                   are the Position of lower left corner and w and h are
%                   the width and height of the axes object eg [.1 .1 .4 .4]
%       Method 2: 	cell array with 3 integers or 3x1 vector used as argument to subplot()
%                   command. e. g. {2 2 [2 4]} or [2 2 2]
% 'Interpreter': interpreter to be used for x/y labels and for title
% 'Title':  plot to specified axes object. (Def. none)
% 'YLabel': the y axes label (Def. none)
% 'XLabel': the x axes label (Def. none)
% 'ZLabel': the z axes label (Def. none)
% 'XLim': the x axes limits (Def. [] meaning limits will be detrmined by matlab)
% 'YLim': the y axes limits of left y axis (Def. [] meaning limits will be detrmined by matlab)
% 'YLimRight': the y axes limits of right y axis (Def. [] meaning limits will be detrmined by matlab)
% 'Hold': Hold axes after plotting. Either on or off (Def. 'on')
% 'Legend': either on or off to display legend (Def. 'on')
% 'Grid': either on or off (Def. 'on')

% Parse Axes related Options
p = inputParser; 
p.KeepUnmatched = true;
addParameter(p,'Position',{1 1 1},@(x) (iscell(x) && isvector(x) && length(x) == 3) || (isnumeric(x) && isvector(x) && (length(x) == 4 || length(x) == 3)));
addParameter(p,'Interpreter','latex',@(x) any(strcmpi(x,{'none' 'tex' 'latex'})));
addParameter(p,'Title','',@ischar);
addParameter(p,'Ylabel','',@ischar);
addParameter(p,'Zlabel','',@ischar);
addParameter(p,'Xlabel','',@ischar);
addParameter(p,'Ylim',[],@(x) isequal(size(x(:)),[2 1]) && x(2)>=x(1) || isempty(x));
addParameter(p,'YlimRight',[],@(x) isequal(size(x(:)),[2 1]) && x(2)>=x(1) || isempty(x));
addParameter(p,'Xlim',[],@(x) isequal(size(x(:)),[2 1]) && x(2)>=x(1) || isempty(x));
addParameter(p,'Hold','off',@(x) any(strcmpi(x,{'on' 'off'})));
addParameter(p,'Legend','on',@(x) any(strcmpi(x,{'on' 'off'})));
addParameter(p,'Grid','on',@(x) any(strcmpi(x,{'on' 'off'})));
parse(p,varargin{:});

axopt = p.Results;
axopt.Unmatched = p.Unmatched;
end


function setAxesOptions(axopt,ax)
%SETAXESOPTIONS Sets Axes Options struct to specified axes object. If ax is
%not supplied then gca will be used.
%
% setAxesOptions(axopt,ax)
% axopt [struct]: Valid axes options structure. See idModels.util.axesOptions.m
% ax [opt. Axes or Axes handle]: Axes object/handle for wchich options will be set.

if nargin == 1
    ax = gca;
end
axes(ax);
if strcmpi(axopt.Legend,'on') 
    legend off; 
    wrg = warning;
    warning('off');
    lgd = legend('show');
    warning(wrg);
    lgd_str = get(lgd,'String');
    set(lgd,'Interpreter',axopt.Interpreter);
    if isempty(lgd_str) || all(cellfun(@isempty,lgd_str))
        legend off;
    end
end
if ~isempty(axopt.Title) title(axopt.Title,'Interpreter',axopt.Interpreter); end
if ~isempty(axopt.Ylabel) ylabel(axopt.Ylabel,'Interpreter',axopt.Interpreter); end
if ~isempty(axopt.Xlabel) xlabel(axopt.Xlabel,'Interpreter',axopt.Interpreter); end
if ~isempty(axopt.Zlabel) zlabel(axopt.Zlabel,'Interpreter',axopt.Interpreter); end
if ~isempty(axopt.Ylim)
    if verLessThan('matlab','9.0') % just use current yaxis since there is no right axes
        ylim(axopt.Ylim);
    else % make sure to set left axis
        ax.YAxis(1).Limits = axopt.Ylim;
    end
end
if ~isempty(axopt.YlimRight) 
    if verLessThan('matlab','9.0')
        warning('Plotting to the right axis is unsupported for matlab versions below 2016a!');
    else
        if length(get(ax,'YAxis'))>1
            ax.YAxis(2).Limits = axopt.YlimRight;
        end
    end
end
if ~isempty(axopt.Xlim) xlim(axopt.Xlim); end
if strcmpi(axopt.Grid,'on') grid on; end
if strcmpi(axopt.Hold,'on') hold on; end
end
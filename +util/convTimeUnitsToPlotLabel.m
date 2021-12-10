function lab = convTimeUnitsToPlotLabel(Tu)
%CONVTIMEUNITSTOPLOTLABEL Maps a Time Unit String to a Plot Label, e. g. 'seconds' is mapped to 's'.
%
% lab = convTimeUnitsToPlotLabel(Tu)
% lab [char]: Plotlabel 
% Tu [TimeUnit]: TimeUnit

switch Tu
    case 'nanoseconds'
        lab = 'ns';
	case 'microseconds'
        lab = '\mu s';
	case 'milliseconds'
        lab = 'ms';
    case 'seconds'
        lab = 's';
    case 'minutes'
        lab = 'min';
    case 'hours'
        lab = 'h';
    case 'days'
        lab = 'd';
    otherwise
        error('Unsupported input!');
end
end

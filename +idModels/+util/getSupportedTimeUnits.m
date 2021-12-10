function allowedUnits = getSupportedTimeUnits()
%Returns Cell array of Strings of supported Timeunits. (Used within TIntervals class).
%
% allowedUnits = getTimeUnits()
% allowedUnits [cell of strings]: 

%% WARNING: If Any changes are made here make sure to ensure compatibility with getTimeUnitConvFactor(timeUnit)
allowedUnits = {'nanoseconds','microseconds','milliseconds','seconds','minutes','hours','days','weeks','months','years'};
end
function fac = getTimeUnitConvFactor(timeUnit)
%Returns the timefactor which allows conversion of a timeUnit (inputargument) into days. 
%
% fac = getTimeUnitConvFactor(timeUnit)
% timeUnit [string]: one of the following: 'nanoseconds','microseconds','milliseconds','seconds','minutes','hours','days','weeks','months','years'
% fac [scalar double]: the conversion factor
% Example: getTimeUnitConvFactor('hours') returns 1/24.


%% WARNING: If Any changes are made here make sure to ensure compatibility with getSupportedTimeUnits()

import time.*
if ~ischar(timeUnit)
   error('The Inputargument "timeUnit" needs to be a String!'); 
end

allowedUnits = idModels.util.getSupportedTimeUnits();
idx_match = strcmpi(allowedUnits,timeUnit);

if ~any(idx_match)
    error('Unsupported Time Unit!');
end
    
factors = [1e-9 1e-6 1e-3 1 60 3600 86400 604800 2629800 31557600]/86400;
fac = factors(strcmpi(allowedUnits,timeUnit));
end


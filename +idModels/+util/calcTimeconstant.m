function T = calcTimeconstant(poles,Ts)
%CALCTIMECONSTANT Transforms an pos. real discrete pole into an timeconstant.
%This function can be used to get a rough idea of speed of converegence 
%of a dicrete time system with given poles. 
%NOTE: Positive real poles are only supported currently. 
%
% T = calcTimeconstant(poles,Ts)
% poles [double vector]: poles of dicrete time system
% Ts [pos. scalar]: Sampling time 
% T [double vector]: associated time constants of ZOH converted continous system.

assert(all(isreal(poles)) && all(poles>=0),'Complex and negative poles arent supported yet!');
T = -Ts./log(poles);
end


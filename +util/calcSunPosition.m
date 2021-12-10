function [elevation, azimuth] = calcSunPosition(time,lambda,phi)
%Berechnet Zeitreihen der Sonneposition in Abhängigkeit der geografischen Lage und der Ortszeit.
%
% [elevation, azimuth] = calcSunPosition(time,lambda,phi)
% Time [datenum Vektor]: !!!CET/UTC+1!!!-Zeit in Matlab-Format. 
% Lambda [double]: Geographische Länge im Gradmaß (W == neg; O == pos.)
% Phi [double]: Geographische Breite im Gradmaß (N == pos; S == neg.)
% Output:
% elevation [vector]: sonnenhöhe in Grad
% azimuth [vector]: azimuth in Grad

%Check inputs
assert(isscalar(lambda) && lambda>=-180 && lambda<=180,'Lambda needs to be a scalar in the range of -180 <= lambda <= 180');
assert(isscalar(phi) && phi>=-90 && phi<=90,'Phi needs to be a scalar in the range of -90 <= phi <= 90');

phi = phi/180*pi; %in Bogenmaß umrechnen

Time = time(:);
t_dv = datevec(Time); % zeit im datevector format
leapyear = mod(t_dv(:,1),4)==0; % indikator ob schaltjahr

% Zeitstempel relativ zum 1.Januar des aktuellen Jahres bestimmen:
DateVec=datevec(Time(1));
Time=(Time-datenum(DateVec(1),1,1))*24;        %Zeit ab hier in h

J=2*pi*Time/24./(365+leapyear);

%Central European time
CET = mod(Time,24);

% Mittlere Ortszeit
MOZ = CET - (15 - lambda)/15;

%Zeitgleichung in Minuten
ZGL = 0.0066 + 7.3525*cos(J + pi/180*(85.9)) + 9.9359*cos(2*J + pi/180*(108.9))...
   + 0.3387*cos(3*J + pi/180*(105.2));

% Wahre Ortszeit
WOZ = MOZ + ZGL./60;

% Stundenwinkel (positiv zum Nachmittag)
omega = (WOZ - 12)*15*pi/180;

% Sonnendeklination in Grad
delta = 0.3948 - 22.2559*cos(J + pi/180*(9.1)) - 0.3915*cos(2*J + pi/180*(5.4))...
   - 0.1764*cos(3*J + pi/180*(26));

% Vektor in Richtung Sonne: mit Komponenten in Richtung [Süden, Osten, Zenit]
RadiationVector = [-cos(phi) .* sin(pi/180*delta) + ...
   cos(pi/180*delta) .* cos(omega) .* sin(phi),-cos(pi/180*delta) .* sin(omega),cos(pi/180*delta) ...
   .* cos(phi) .* cos(omega) + sin(pi/180*delta) .* sin(phi)];

% Sonnenelevation
gamma = asin(cos(omega) .* cos(phi) .* cos(pi/180*delta) + sin(phi) .* sin(pi / 180*delta));

elevation = gamma*180/pi;
azimuth=(pi+atan2(-RadiationVector(:,2),RadiationVector(:,1)))*180/pi;

end

function [Edir_gen, Ediff_gen, Eref_gen] = calcRadiationOnInclinedPlane(Edir,azimuth,height,alpha_gen,gamma_gen,Ediff,method)
%Berechnet die Strahlung auf einer geneigten Ebene welche um einen azimuthwinkel alpha_gen gegenüber 
%Norden gedreht und einen Höhenwinkel gamma_gen gegenüber der horizontalen Erdoberfläche angekippt ist und dem Sonnenstandsvektor
%dessen  Position durch azimuth und höhe beschrieben wird. Die routine berücksichtigt den Fall für welchen der Winkel zw. sonne
%und Normalen der Ebene >90 Grad werden. Dann steht die Sonne hinter der ebene und die Strahlung wird daher zu 0 gesetzt.
%Edir_gen ist stets größer oder gleich Edir.  
%
% [Edir_gen, Ediff_gen] = calcEffectiveRadiationOnInclinedPlane(Edir,azimuth,height,alpha_gen,gamma_gen,Ediff)
% elevation [vector]: Sonnenhöhe in Grad
% azimuth [vector]: Azimuth der Sonne in Grad
% alpha_gen [pos. double]: Azimuthwinkel in Grad um welche die Ebene gegüber Norden gedreht ist. 270 würde Westausrichtung entsprechen. 90 Ost.
% gamma_gen [double]: Höhenwinkel gamma_gen gegenüber der horizontalen Erdoberfläche in grad. 0 entspricht einer horiz. Ebene.
% Edir [vector]: Zeitreihe der Direktstrahlung auf der Horizontalen in W/m^2
% Ediff [opt. vector]: Zeitreihe der Diffusstrahlung auf der Horizontalen in W/m^2
% method [int 1 - 3]: Methode um den Diffusstrahlungsanteil auf der gen. Ebene zu berechnen.
%   Siehe V. Quasching: "Solare Energiesysteme" (2013), Abschn. 2.6.2 
%   1: isotroper Ansatz (nimmt Ggleichverteilte Strahlungsdichte an). Relativ ungenau bei bedecktem Himmel
%   2: Klucher- Modell 1979
%   3: Perez Model

%assert(all(size(Edir)==size(azimuth)) && all(size(Edir)==size(height)))
height(height<5) = 0; % Die Umrechnung wird sehr ungenau wenn die Sonnenstrahlen parallel zum Horizont

if nargin <= 6
    method = 2;
end

% In Bogenmaß umrechnen
gamma_gen = gamma_gen/180*pi;
alpha_gen = alpha_gen/180*pi;
az = azimuth/180*pi;
h = height/180*pi;

% Für kart. KS: x: Norden y: Westen z: Zenitrichtung gilt für den Vektor der Sonne s und den Normalenvektor der Ebene n
% s = [cos(h)cos(az); -cos(h)sin(az); sin(h)];
% n = [sin(gamma_gen)cos(alpha_gen); -sin(gamma_gen)sin(alpha_gen); cos(gamma_gen)];
% damit ergibt sich für cosTheta = sn (Additionstheorem beachten)
cosTheta = cos(h).*sin(gamma_gen).*cos(az-alpha_gen) + sin(h)*cos(gamma_gen); % cos des winkels zw. Sonnenvektor und Normalenvektor der gen Ebene

Edir_gen = Edir.*cosTheta./sin(h); % Strahlung auf gen. Ebene
Edir_gen(cosTheta<0 | h <= 0) = 0; % wenn Sonne hinter Ebene oder sonne untergegangen

if nargout > 1 
    % verschiedene Modelle: siehe Volker Quaschning - Regenerative Energiesysteme S.73 ff.
    if method == 1
        Ediff_gen = 0.5*Ediff*(1+cos(gamma_gen));
    elseif method == 2
        Eglob = Ediff+Edir;
        F = 1-(Ediff./Eglob).^2;
        F(Eglob==0) = 0;
        Ediff_gen = 0.5*Ediff.*(1+cos(gamma_gen)).*(1+F*sin(gamma_gen/2)^3).*(1+F.*cosTheta.^2.*cos(h).^3);
    elseif method == 3
        kappa = 1.041;
        theta_zen = pi/2 - h; % solarer zenith winkel
        AM = 1./cos(theta_zen); 
        AM(AM>30) = 30; %https://pvpmc.sandia.gov/modeling-steps/1-weather-design-inputs/irradiance-and-insolation-2/air-mass/
        E0 = 1360.8; %[E/m^2] Solarkonstante 
		Delta = AM.*Ediff/E0;
        F = [-0.008	0.130   0.330   0.568   0.873   1.132   1.060   0.678
            0.588 	0.683   0.487   0.187   -0.392  -1.237  -1.6    -0.327
            -0.062	-0.151  -0.221  -0.295  -0.362  -0.412  -0.359  -0.25
            -0.06   -0.019  0.055   0.109   0.226   0.288   0.264   0.156
            0.072   0.066   -0.064  -0.152 	-0.462  -0.823  -1.127 -1.377
            -0.022  -0.029 	-0.026  -0.014  0.001   0.056   0.131   0.251];
        epsilon = ((Ediff + Edir./sin(h))./Ediff + kappa*theta_zen.^3)./(1 + kappa*theta_zen.^3);
        epsilon(isnan(epsilon)) = 8; % 0/0
        Bin = 8*ones(length(epsilon),1); 
        Bin(epsilon>=1 & epsilon<1.065) = 1;
        Bin(epsilon>=1.065 & epsilon<1.23) = 2;
        Bin(epsilon>=1.23 & epsilon<1.5) = 3;
        Bin(epsilon>=1.5 & epsilon<1.95) = 4;
        Bin(epsilon>=1.95 & epsilon<2.8) = 5;
        Bin(epsilon>=2.8 & epsilon<4.5) = 6;
        Bin(epsilon>=4.5 & epsilon<6.2) = 7;
        F1 = F(1,Bin)' + (F(2,Bin)').*Delta + (F(3,Bin)').*theta_zen;
        F2 = F(4,Bin)' + (F(5,Bin)').*Delta + (F(6,Bin)').*theta_zen;
        a = max(0,cosTheta);
        b = max(0.0872,cos(theta_zen));
        Ediff_gen = Ediff.*( (1-F1).*(1+cos(gamma_gen)/2 + a.*F1./b + F2.*sin(gamma_gen)));
        Ediff_gen(Ediff_gen<0) = 0;
    end
end

if nargout > 2 
    A = 0.2; % Albedo wert
    Eref_gen = (Edir + Ediff)*A*0.5*(1-cos(gamma_gen)); 
end
end

%!test
%! az_s = [270 90 270]; h_s = [45 45 90]; az_w = 270; h_w = 45; 
%! Edir_gen = util.calcRadiationOnInclinedPlane(ones(1,length(az_s)),az_s,h_s,az_w,h_w,1);
%! assert(abs(Edir_gen(1) - 1/cos(45/180*pi)) <= eps && Edir_gen(2)==0 && Edir_gen(3)==cos(45/180*pi))

% %! rad = 'fhg_eas...drewag.o..global.radiation';
% %! phi = 51 + 2/60;
% %! lambda = 13 + 44/60;
% %! tsc = signals.TSignalCollection.loadFromDB('url','jdbc:postgresql://vm020:5431/raum51','user','eas_domain','pw','masd.ZghA_23',...
% %!                                       'timespan',{'2016-4-1' '2016-9-1'},'identifiers',rad);
% %! [hd, azd] = util.calcSunPosition(tsc(rad).getabstime(1),lambda,phi)
% %! tsc.sun_height = tsc(rad); tsc.azimuth = tsc(rad);
% %! tsc.sun_height.Data = hd; tsc.azimuth.Data = azd;
% %! tsc.Edir_gen = tsc(rad).calcRadiationOnInclinedPlane(tsc.azimuth,tsc.sun_height,az_w,h_w);
% %! tsc.plot('Identifiers',{{'Edir_gen' rad} {'sun_height' 'azimuth'}})
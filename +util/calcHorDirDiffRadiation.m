function [Edir,Ediff] = calcHorDirDiffRadiation(EglobHor,Hsun,varargin)
%CALCHORDIRDIFFRADIATION Seperates Global horizontal radiation to diffuse and direct part using
%model proposed by Reindl, D. T. and Beckman, W. A. and Duffie, J. A. "Diffuse Fraction 
%Correlations" (1990). Alternatively see also Quasching, V. "Solare Energiesysteme" (2013) p. 67.
%
% [Edir,Ediff] = calcHorDirDiffRadiation(EglobHor,Hsun)
% Edir [double vector]: Direct Radiation on horizontal plane [W/m2]
% Ediff [double vector]: Diffuse Radiation on horizontal plane [W/m2]
% EglobHor [double vector]: Global Radiation on horizontal plane [W/m2]
% Hsun [double vector]: Height of sun [degress]

%Hsun(Hsun<1) = 1;
Hs = Hsun*pi/180;
sinHsun = sin(Hs);
E0 = 1360.8; %[E/m^2] Solarkonstante 
kT = EglobHor./(E0*sinHsun);
Ediff = EglobHor.*(1.4 - 1.749*kT + 0.177*sinHsun);
ix = kT<=0.3;
Ediff(ix) = EglobHor(ix).*(1.02 - 0.254*kT(ix) + 0.0123*sinHsun(ix));
ix = kT>=0.78;
Ediff(ix) = EglobHor(ix).*(0.486*kT(ix) - 0.182*sinHsun(ix));
ix = Ediff>EglobHor;
Ediff(ix) = EglobHor(ix);
Edir = EglobHor - Ediff;
%figure; plot([EglobHor Edir Ediff]); legend({'Global' 'Direkt' 'Diffus'})
end


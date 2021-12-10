function [f,df] = fun_rad(u,p)
%Calculates effective radiation through the windows. Three methods are available: 
% length(p) == 2:   uses global horizontal radiation and converts it to 
%                   global rad on inclined plane
% length(p) == 3:   Splits global rad to diffuse and direct part and
%                   converts to dir, diff and diff reflected radiation on 
%                   inclined plane using method specified by p(3) (fixed parameter!).
% length(p) == 5:   Like length(p) == 3 but direct part is multiplied by
%                   logistic function. p(4:5) are used to scale logistic
%                   function.
%u(1,:) -> global rad on hor plane [!!!W/m^2!!!] (do not use)
%u(2,:) -> azimuth [deg]
%u(3,:) -> sun height [deg]
%u(4,:) -> shading position (if size(u,2)>3) 0 ... opened, 1 ... closed
%p(1) -> rotation angle from north (west->270) 
%p(2) -> horizontal orientation (perpendicular to surface of earth ->90)
%p(3) -> (opt.) method for computing diffuse radiation on inclined plane. see util.calcRadiationOnInclinedPlane
%p(4:5) -> parameters used to scale and shift logistic function. see idModels.func.fun_logistic

if length(p) == 2 
    f = util.calcRadiationOnInclinedPlane(u(:,1),u(:,2),u(:,3),p(1),p(2));
    if nargout > 1 
        df = NaN(length(f),2); %Keine Lust Gradienten zu berechnen... Bisher auch noch nicht nötig, da keine Optimierungsparameter
    end
else
    [Edir,Ediff] = util.calcHorDirDiffRadiation(u(:,1),u(:,3));
    [Edir_gen,Ediff_gen,Erefl_gen] =  util.calcRadiationOnInclinedPlane(Edir,u(:,2),u(:,3),p(1),p(2),Ediff,p(3));
    df = NaN(size(u,1),3); % Keine Lust Gradienten zu berechnen...
    if length(p) == 5
        [L,dL] = idModels.func.fun_logistic(u(:,3),p(4:5));
        Edir_gen = L.*Edir_gen; 
        df = [df dL.*Edir_gen];
    end
   	f = Edir_gen + Ediff_gen + Erefl_gen;
end
if size(u,2)>3
    f = (1-u(:,4)).*f;
 	if nargout > 1 && length(p)>1
        df = (1-u(:,4)).*df;
  	end
end
end
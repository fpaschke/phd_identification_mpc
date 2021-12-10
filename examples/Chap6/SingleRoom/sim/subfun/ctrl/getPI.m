function [Kp,Ki] = getPI(Gin,varargin)

p = inputParser();
addParameter(p,'DiscMth','zoh');
addParameter(p,'Setting','Fast',@(x) any(strcmpi(x,{'fast' 'slow' 'aper_dist' 'aper_ref' '20_dist' '20_ref'})));
addParameter(p,'PlotStepResponse',false,@(x) x == 1 || x == 0)
addParameter(p,'Method','Tsum',@(x) any(strcmpi(x,{'chr' 'tsum' 'absopt'}))); 
addParameter(p,'Ts',Gin.Ts)
parse(p,varargin{:});

P = Gin;
assert(isequal(Gin.TimeUnit,'seconds'))
if P.Ts>0
   P = d2c(P,p.Results.DiscMth); 
end
if ~isa(P,'zpk')
   P = zpk(P); 
end
Ks = P.K*prod(-P.Z{1})/prod(-P.P{1}); 
T = -1./P.P{1};     % TimeConstants in Gin.TimeUnit
switch lower(p.Results.Method)
    case 'tsum'
        Td = -1./P.Z{1};    % Differentiating TimeConstants in Gin.TimeUnit
        Tsum = sum(T) - sum(Td); % 
        % R(s) = Kr*(1+1/Tn*s) -> Ki = Kr/Tn
        if strcmpi(p.Results.Setting,'fast')
            Kr = 1/Ks;
            Tn = 0.7*Tsum;
        elseif strcmpi(p.Results.Setting,'slow')
            Kr = 0.5/Ks;
            Tn = 0.5*Tsum;   
        end
    case 'chr'
        t = 0:Gin.Ts:24*3600*20;
        [y,t] = step(P,t);
        dy = [NaN; diff(y)./diff(t)]; 
        [~,ixm] = max(dy);
        m = dy(ixm); n = y(ixm) - m*t(ixm);
        Tu = -n/m; Tg = (Ks-n)/m;
        %plot(t,y); hold on; plot(t(1:150),m*t(1:150)+n,'-r');
        switch lower(p.Results.Setting)
            case 'aper_dist' 
                Kr = 0.6*Tg/(Tu*Ks); Tn = 4*Tu;
            case 'aper_ref' 
                Kr = 0.35*Tg/(Tu*Ks); Tn = 1.2*Tg;
            case '20_dist' 
                Kr = 0.7*Tg/(Tu*Ks); Tn = 2.3*Tu;
            case '20_ref'
                Kr = 0.6*Tg/(Tu*Ks); Tn = 1*Tg;
            otherwise % aper ref
                Kr = 0.35*Tg/(Tu*Ks); Tn = 1.2*Tg;
        end 
    case 'absopt'
        assert(length(P.P{1})==2 && all(isreal(P.P{1})),'Only implemented for 2nd order systems!');
        Tn = max(T);
        Kr = Tn/(2*Ks*min(T));
end
% METHOD 1: Discretization of integrator using Backward Euler 1/s = Ts/(1-z^-1)
% Ki = Kr/(Tn)*Gin.Ts; %1/(K*s)*s = 1/K 
% Rd = tf(Kr,1,Gin.Ts) + tf([Kr/Tn 0],[1 -1],Gin.Ts); Rd.Variable = 'z^-1';

% METHOD 2: Use c2d() to get R(z) = b0 + b1 z^-1/(1-z^-1) and recover cobntroller parameters
Rd = c2d(tf([Kr*Tn Kr],[Tn 0]),p.Results.Ts,p.Results.DiscMth); Rd.Variable = 'z^-1';
% Intended controller: R = Kr + Ki*1/1-z^-1 = (Kr+Ki - Kr z^-1)/(1-z^-1) 
Kp = -Rd.Numerator{1}(2);
Ki = Rd.Numerator{1}(1) + Rd.Numerator{1}(2); 
G = minreal(zpk(Rd*Gin/(1+Rd*Gin)));
if ~isstable(G)
    warning('Closed loop isnt stable!');
end
if p.Results.PlotStepResponse
    G = zpk(minreal(G,1e-8));
    step(chgTimeUnit(G,'hours'));
end
end


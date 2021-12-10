function [Tlower, Tupper] = buildTsetSig(t,SigOpt,Tout)
%BUILDTSETSIG Build signal for upper and lower Temperature bounds.

isnoday = rem(t,7*86400) < 2*86400 | rem(t,86400) < SigOpt.DayTime(1)*3600 | rem(t,86400) >= SigOpt.DayTime(2)*3600;
N = length(t);
if ~isfield(SigOpt,'T_set_upper')
    SigOpt.T_set_upper =  SigOpt.T_set_lower;
end
Nz = size(SigOpt.T_set_lower,1);

if nargin > 2
    NsDay = 86400/diff(t(1:2));
    Nday = floor(N/NsDay);
    Tout = reshape(Tout(1:NsDay*Nday),NsDay,Nday);
    ToutMax = kron(max(Tout,[],1)',ones(NsDay,1));
end

for nz = 1:Nz
    Tlower(:,nz) = SigOpt.T_set_lower(nz,2)*ones(N,1); 
    Tlower(isnoday,nz) = SigOpt.T_set_lower(nz,1); 
  	Tupper(:,nz) = SigOpt.T_set_upper(nz,2)*ones(N,1); 
    Tupper(isnoday,nz) = SigOpt.T_set_upper(nz,1); 
    if nargin > 2
       	Tup = interp1([-100; SigOpt.T_out_summer(:); 100]',[SigOpt.T_set_upper(nz,2) SigOpt.T_set_upper(nz,2) SigOpt.T_set_summer SigOpt.T_set_summer]',ToutMax);
        Tup(NsDay*Nday+1:N) = Tup(end);
        Tup = round(2*Tup)/2;
        Tupper(~isnoday,nz) = Tup(~isnoday); 
    end
    if isfield(SigOpt,'PrbsAmplitude')        
        rng(nz);     % Seed for reproducebility
        Nw = (t(end) - t(1))/(7*86400); 
        Ts = diff(t(1:2))/60;
        Trnd = circshift(kron(2*SigOpt.PrbsAmplitude(nz)*(round(rand(ceil(Nw*7*24/SigOpt.PrbsInterval(nz)),1))-0.5),ones(SigOpt.PrbsInterval(nz)*60/Ts,1)),...
            (nz-1)*(SigOpt.PrbsPhaseShift*60/Ts));
     	if isfield(SigOpt,'PrbsDayOnly') && SigOpt.PrbsDayOnly
            ixd = repmat(   [false(2*24*60/Ts,1)
                            repmat([false(floor(SigOpt.DayTime(1)*60/Ts),1)
                                    true(ceil(diff(SigOpt.DayTime)*60/Ts),1) 
                                    false((24-SigOpt.DayTime(2))*60/Ts,1)],5,1)],Nw,1);
            Trnd = Trnd.*ixd;
        end
        Trnd = [Trnd; zeros(length(t)-length(Trnd),1)];
        Tlower(:,nz) = Tlower(:,nz) + Trnd(1:size(Tlower,1));
        Tupper(:,nz) = Tupper(:,nz) + Trnd(1:size(Tlower,1));
    end
end
end
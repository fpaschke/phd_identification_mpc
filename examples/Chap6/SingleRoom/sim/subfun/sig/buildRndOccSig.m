function O = buildRndOccSig(t,Max,Schedule)
%BUILDOCCSIG Builds random occupancy signal.

N_w = (t(end)-t(1))/(86400*7);
T_s = unique(diff(t))/60;
Ns = diff([0 Schedule 24])*60/T_s;
O = [];
for nw = 1:N_w
    for d = 1:7
       if d == 1 || d == 2
           O = [O; zeros(1440/T_s,1)];
       else
           for ns = 1:length(Ns)
               if ns == 1 || ns == length(Ns)
                 	O = [O; zeros(Ns(ns),1)];
               else
                    O = [O; round(Max*rand)*ones(Ns(ns),1)];
               end
           end
       end
    end
end
wd = weekday(t(1)/86400);
O = circshift(O,(6-wd)*(1440/T_s));
O = [O; O(end)];
end


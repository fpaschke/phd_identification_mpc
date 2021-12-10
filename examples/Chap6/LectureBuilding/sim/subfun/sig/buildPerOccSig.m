function O = buildPerOccSig(Nw,Ts,M,S,T,Tr,Seed)
%BUILDOCCSIG Build an Occupancy Signal.
%
% function O = buildOccSig2(Nw,Ts,Nmax,Std,T,Tr,Seed)
% Nw [pos int]:             Number of weeks 
% Ts [pos int]:             Sampling Time in min
% M [double 5 x nh]:      	Mean number of people. Row corresponds to day, column to number of hour
% S [double]:               Standard dev. of number of people.
% T [opt. N x 2 double]:	Start end End-Time of Occupancy at day in hours (schedule).
%                         	Def. T = [7.5 9; 9.5 11; 11.5 13; 13.5 15;15.5 17; 17.5 19];
% Tr [opt double]:          Rise/Falling Time in min
% Seed [scalar double]:   	Seed for Random Number generator

if nargin <7 || isempty(Seed); Seed = 0; end
rng(Seed); % For Reproducability                  	
if nargin<5 || isempty(T)
    T = [   7.5 9; 
            9.5 11; 
            11.5 13;
            13.5 15;
            15.5 17;
            17.5 19];
end
if nargin<6 || isempty(Tr) || Tr == 0
    Tr = 1e-8;
end
Nh = size(T,1);
if size(M,1)==1; M = M.*ones(5,1); end
if size(M,2)==1; M = M.*ones(1,Nh); end
assert(all(diff(T,1,2)>0) && all(T(:)<=24));

% Build Daysignal
T = 60*T; % In Min
T_break = [0; diff([T(1:end-1,2) T(2:end,1)],[],2)];
t = []; v = [];
for nd = 1:7
    t = [t (nd-1)*1440]; 
    v = [v 0];
    if nd > 2
        for nh = 1:Nh
            if t(end)>=(T(nh,1)-Tr)+(nd-1)*1440
                t(end) = mean([T(nh-1,2) T(nh,1)]) + (nd-1)*1440;
                v(end) = M(nd-2,nh)*(1 - T_break(nh)/(2*Tr));
                t = [t [T(nh,1) T(nh,2) T(nh,2)+Tr] + (nd-1)*1440]; 
                v = [v M(nd-2,nh) M(nd-2,nh) 0];
            else
                t = [t [T(nh,1)-Tr T(nh,1) T(nh,2) T(nh,2)+Tr] + (nd-1)*1440]; 
                v = [v 0 M(nd-2,nh) M(nd-2,nh) 0];
            end
        end
    end
end
t = [t 1440*7-1]; v = [v 0];
w = interp1(t,v,0:Ts:1440*7-1)';
O = repmat(w,Nw,1);
O = [O; 0];
O = round(O + O.*randn(length(O),1)*S);
O(O<0) = 0;
end
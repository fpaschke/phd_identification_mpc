function resample(obj,Ts,mth)
%RESAMPLE Resamples Model to a new Sampling time Ts using specified Method.
%
% obj = resample(obj,Ts,method)
% obj [NsfPolyModel]: The Model to be Resampled.
% Ts [pos. double scalar]: The new Sampling Time in obj.TimeUnit
% method [opt. char]: either 'zoh' or 'tustin' 

if nargin < 3 || isempty(mth) 
    mth = 'zoh'; 
end % to preserve stability
assert(isscalar(Ts) && Ts >0,'Ts need to be a pos. scalar!');  
warning off;
obj.defactorize();
warning on;
oldInfo = obj.Info;


%% Using Ctrl Toolbox
try
	assert(isdiag(obj.NoiseVariance),'"NoiseVarance" isnt a diag Matrix! -> Unsupported yet!');
    [G_num,G_den,H_num,H_den] = obj.calcGH('AbsorbNoiseVariance',false);
  	sys = d2d(tf([G_num H_num],[G_den H_den],obj.Ts,'Variable','z^-1'),Ts,mth); % resample system   
    nu = size(G_num,2);
    %F = calcKStepPredictor(obj,ceil(Ts/obj.Ts));
    for ny = 1:obj.OutputDimension
        % find poles which are common to all Gi and H
        ix_num_0 = cellfun(@(n) all(n==0), sys.num(ny,:)); % indexes of tf where numerator is 0
        Rem = repmat({1},1,size(sys,2));
        [Ai,Rem(~ix_num_0)] = obj.commonZeros(sys.den(ny,~ix_num_0),true);
        newA(ny,ny).val = Ai;
        newA(ny,ny).free = [false true(1,length(Ai)-1)];
        newA(ny,ny).std = NaN(1,length(newA(ny,ny).val));
        newA(ny,ny).factorized = false;
        
        % find H
      	ixn = size(sys,2)-obj.OutputDimension+ny;
        newD(ny,1).val = Rem{ixn};
        newD(ny,1).free = [false true(1,length(newD(ny,1).val)-1)];
        newD(ny,1).std = NaN(1,length(newD(ny,1).val));
        newD(ny,1).factorized = false;
        
        % find Poles which are common to all Gi
        [Ei,Rem(~ix_num_0(1:nu))] = obj.commonZeros(Rem(~ix_num_0(1:nu)),true);
      	newE(ny,1).val = Ei;
        newE(ny,1).free = [false true(1,length(Ei)-1)];
        newE(ny,1).std = NaN(1,length(newE(ny,1).val));
        newE(ny,1).factorized = false;
        
        % build B and F
        for ni = 1:nu
            newB(ny,ni).val = sys.num{ny,ni};
            newB(ny,ni).free = true(1,length(newB(ny,ni).val));
            newB(ny,ni).std = NaN(1,length(newB(ny,ni).val));
            newB(ny,ni).factorized = false;
            newF(ny,ni).val = Rem{ni};
            newF(ny,ni).free = [false true(1,length(newF(ny,ni).val)-1)];
            newF(ny,ni).std = NaN(1,length(newF(ny,ni).val));
            newF(ny,ni).factorized = false;
        end
        
        % Noise
        %obj.NoiseVariance(ny,ny) = sys.num{ny,ixn}(1)^2;
        newC(ny,1).val = sys.num{ny,ixn}./(sys.num{ny,ixn}(1));
        newC(ny,1).free = [false true(1,length(newC(ny,1).val)-1)];
       	newC(ny,1).std = NaN(1,length(newC(ny,1).val));
        newC(ny,1).factorized = false;
        newLambda(ny,ny) = obj.NoiseVariance(ny,ny)*sys.num{ny,ixn}(1);%*sum(F{ny,ny}.^2);
    end
%% Using SysId Toolbox
catch
    A = arrayfun(@(Ai) Ai.val,obj.A,'UniformOutput',0);
    B = arrayfun(@(Bi) Bi.val,obj.B,'UniformOutput',0);
    C = arrayfun(@(Ci) Ci.val,obj.C,'UniformOutput',0);
    D = arrayfun(@(Di) Di.val,obj.D,'UniformOutput',0);
    F = cell(obj.OutputDimension,size(obj.B,2));
    for no = 1:obj.OutputDimension
       	C{no} = conv(C{no},obj.PreFilter(no).den);
        D{no} = conv(D{no},obj.PreFilter(no).num);
        F(no,:) = arrayfun(@(Fi) conv(Fi.val,obj.E(no).val),obj.F(no,:),'UniformOutput',0);
    end
    Sys = d2d(idpoly(A,B,C,D,F,obj.NoiseVariance,obj.Ts),Ts,mth);
    for ny = 1:obj.OutputDimension
        for ny2 = 1:obj.OutputDimension
            if iscell(Sys.A)
                newA(ny,ny2).val = Sys.A{ny,ny2};
                newC(ny,1).val = Sys.C{ny};
                newD(ny,1).val = Sys.D{ny};
            else
                newA(ny,ny2).val = Sys.A;
                newC(ny,1).val = Sys.C;
                newD(ny,1).val = Sys.D;
            end
            newA(ny,ny2).free = [false true(1,length(newA(ny,ny2).val)-1)];
            newA(ny,ny2).std = NaN(1,length(newA(ny,ny2).free));
         	newC(ny,1).free = [false true(1,length(newC(ny).val)-1)];
            newC(ny,1).std = NaN(1,length(newC(ny,1).free));
            newD(ny,1).free = [false true(1,length(newD(ny).val)-1)];    
            newD(ny,1).std = NaN(1,length(newD(ny,1).free));
        end
        for ni = 1:size(B,2)
            if iscell(Sys.B)
                newB(ny,ni).val = Sys.B{ny,ni};
                newF(ny,ni).val = Sys.F{ny,ni};
            else
                newB(ny,ni).val = Sys.B;
                newF(ny,ni).val = Sys.F;
            end
            newB(ny,ni).free = [true(1,length(newB(ny,ni).val))];
            newB(ny,ni).std = NaN(1,length(newB(ny,ni).free));
            newF(ny,ni).free = [false true(1,length(newF(ny,ni).val)-1)];
            newF(ny,ni).std = NaN(1,length(newF(ny,ni).free));
        end
        newE(ny,1).val = 1; 
        newE(ny,1).free = [false];
        newE(ny,1).std = NaN;
    end
    newLambda = Sys.NoiseVariance;
end
%% Set new Polynomials and Ts
obj.Ts = Ts; obj.NoiseVariance = newLambda; obj.Info = [];
obj.A = newA; obj.B = newB; obj.C = newC; obj.D = newD; obj.E = newE; obj.F = newF;
end
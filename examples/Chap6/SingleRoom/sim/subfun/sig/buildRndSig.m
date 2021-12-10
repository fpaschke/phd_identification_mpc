function O = buildRndSig(N_w,T_s,P,Seed)
%BUILD Random input signal.

rng(Seed)
O = zeros(P.TimeSpan(1)*1440/T_s,1)+P.Range(1);
Ns = ceil(diff(P.TimeSpan)*24/P.SwitchInterval);
R = round((diff(P.Range)*rand(Ns,1)+P.Range(1))/P.StepSize)*P.StepSize;
O = [O;
     kron(R,ones(P.SwitchInterval*60/T_s,1))];
Ns = N_w*7*1440/T_s+1;
O = [O;
     zeros(Ns-length(O),1)+P.Range(1)];
if length(O)>Ns
    O = O(1:Ns);
end
end


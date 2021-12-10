function initParameters(M1,M2,varargin)
%Sets all Parametervalues of M1 from M2. This is only possible if all polynomial orders of M2 are less or equal then these 
%of M2 (M1.Na >= M2.Na, M1.Nb >= M2.Nb ... )
%
% initParameters(M1,M2)
% M1 [NsfPolyModel]: Model which parameters will be set
% M2 [NsfPolyModel]: Model from which paramerters will be set.

assert(isequal(size(M1.A),size(M2.A)) && isequal(size(M1.B),size(M2.B)),'The Models are not compatible since their number of Inputs or Outputs do not match!'); 
assert(all(M1.Na(:)>=M2.Na(:)) && all(M1.Nb(:)>=M2.Nb(:)) && all(M1.Nc>=M2.Nc) && all(M1.Nd>=M2.Nd) ...
            && all(M1.Ne>=M2.Ne) && all(M1.Nf(:)>=M2.Nf(:)),'All polynomial orders of M1 need to >= polynomial orders of M2!'); 
assert(all(M2.Nk(:)>=M1.Nk(:)),'Number of Delays Nk in M2 needs to be >= Nk in M1!');

locSetPolyVals('A');
locSetPolyVals('B');
locSetPolyVals('C');
locSetPolyVals('D');
locSetPolyVals('E');
locSetPolyVals('F');
M1.NoiseVariance = M2.NoiseVariance; 
for ni = 1:size(M1.B,2)
    assert(isequal(M1.InputNonlinearity(ni).fun,M2.InputNonlinearity(ni).fun),'Inputnonlinearity needs to be equal!');
    M1.InputNonlinearity(ni).parameters = M2.InputNonlinearity(ni).parameters;
end

function P = locSetPolyVals(P)
    [Nr,Nc] = size(M1.(P));
    for nr = 1:Nr
        for nc = 1:Nc
            hasChanged = locSetEqualParameterization(P,[nr nc]);
            %M1.A(nr,nc).val(2:end) = M*rand(1,length(M1.A(nr,nc).val)-1);
            %M1.A(nr,nc).val(1:length(M2.A(nr,nc).val)) = M2.A(nr,nc).val; 
            M1.(P)(nr,nc).val(1:length(M2.(P)(nr,nc).val)) = M2.(P)(nr,nc).val;
            M1.(P)(nr,nc).val(length(M2.(P)(nr,nc).val)+1:end) = 0;
            locResetParameterization(P,[nr nc],hasChanged)
        end
    end
end

function hasChanged = locSetEqualParameterization(P,ix)
    hasChanged = false;
    if M1.(P)(ix(1),ix(2)).factorized && ~M2.(P)(ix(1),ix(2)).factorized
		M2.factorize(P,ix); 
        hasChanged = true;
    elseif ~M1.(P)(ix(1),ix(2)).factorized && M2.(P)(ix(1),ix(2)).factorized
        M2.defactorize(P,ix); 
        hasChanged = true;
    end
end

function locResetParameterization(P,ix,changeParameterization)
    if changeParameterization
        if ~M2.(P)(ix(1),ix(2)).factorized
            M2.factorize(P,ix); 
        else
            M2.defactorize(P,ix); 
        end
    end
end
end
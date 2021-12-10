function defactorize(obj,polystr,ix)
%DEFACTORIZE Converts polynomial structure P from factorized parameterization to standard form, i. e.:
%P(q) =  K*(q-^1 - p10) * ... * (q-^1 - pn0) = p0 + p1 q^-1 + ... + pn q^-n
%The parametervector val of polynomial P is defined as 
%Factorized description: [K p10 ... pn0].
%Standard (defactorized) description: [p0 p1 ... pn].
%
% defactorize(obj,polystr,ix)
% obj [NsfPolyModel]:   The Model object
% polystr [opt. char]:  	Character that defines the polynomial that should be defactorized ('A' 'B' ... 'F'). If no argument 
%                       is supplied then all polynomials will be defactorized
% ix [opt. 2 x 1 pos. int]: Element of poly that should factorized (polystr need to be supplied).
%
% EXAMPLES:
%
% M = idModels.NsfPolyModel(2,2,2,2);
% Afree = [0 0 0]; A = conv([1 -.9],[1 -.8]);
% M.A.val = A; M.A.free = Afree;
% Bfree = [0 0 1 1]; B = [0 0 5 3]; 
% M.B.val = B; M.B.free = Bfree;
%
% Example 1: Factorize Polynomials A and B and defactorize them again:
% M.factorize('A',[1 1]); M.defactorize('A',[1 1]);
% assert(all(Afree == M.A.free) && all(A == M.A.val))
% M.factorize('B'); M.defactorize('B');
% assert(all(Bfree == M.B.free) && all(B == M.B.val))
%
% Example 2: Factorize all Polynomials and defactorize them again:
% Cfree = M.C.free; C = M.C.val;
% M.factorize()
% M.defactorize()
% assert(all(Afree == M.A.free) && all(A == M.A.val))
% assert(all(Bfree == M.B.free) && all(B == M.B.val))
% assert(all(Cfree == M.C.free) && all(isequaln(C,M.C.val)))

uf = obj.UpdateFlag;
obj.UpdateFlag = false;
if nargin==1 || isempty(polystr)
    obj.defactorize({'A' 'B' 'C' 'D' 'E' 'F'});
    return;
end

if iscell(polystr)
    for nP = 1:length(polystr)
        obj.defactorize(polystr{nP});
    end
    return;
end

P = obj.(polystr);
nk = obj.getPolyDelay(P);
if nargin<3
    [nr,nc] = size(P);
    rows = 1:nr;
    cols = 1:nc;
else
    rows = ix(1);
    cols = ix(2);
end

for r = rows
    for c = cols
        if P(r,c).factorized
            if length(P(r,c).val)>21
                warning('The polynomial %s(%i,%i) is of high order (>20)! Results may become inaccurate!',polystr,r,c);
            end
            P(r,c).val = obj.getCoeff(P(r,c));
            assert(all(imag(P(r,c).val)<obj.ImagTol),'Defactorization resulted in complex valued coefficients!');
            if all(P(r,c).free(nk(r,c)+2:end)) %% all poles are free
                pfree = true(1,length(P(r,c).val)-1-nk(r,c));
            elseif all(~P(r,c).free(nk(r,c)+2:end)) %% all poles are fixed
                pfree = false(1,length(P(r,c).val)-1-nk(r,c));
            else
                warning('Unique transformation of %s(%i,%i).free is not possible. Non-delay parameters will be free after defactorization!',polystr,r,c);
                pfree = true(1,length(P(r,c).val)-1-nk(r,c));
            end
            Pf0 = P(r,c).free(1);
            P(r,c).free = [false(1,nk(r,c)) pfree]; %nk fixed delays, ?, 
            if length(P(r,c).free) ~= length(P(r,c).val)
                P(r,c).free = [Pf0 P(r,c).free];
            end
            P(r,c).std = NaN(1,length(P(r,c).val));
            P(r,c).factorized = false;
        else
            % warning('The Polynomial %s(%i,%i) is in standard form already! -> Skipping!',polystr,r,c);
        end
    end
end
obj.DoChecks = false; 
obj.(polystr) = P; % set factorized status (otherwise checkPoly will fail.)
obj.DoChecks = true;
obj.(polystr) = P; % Perform Checks
obj.UpdateFlag = uf;
if isfield(obj.Info,'Popt'); obj.Info = rmfield(obj.Info,'Popt'); end
if isfield(obj.Info,'CovP'); obj.Info = rmfield(obj.Info,'CovP'); end
end
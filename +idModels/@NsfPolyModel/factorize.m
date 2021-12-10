function factorize(obj,polystr,ix)
%FACTORIZE Converts polynomial structure P(q^-1) to factorized parameterization, i. e.:
%P(q) = p0 + p1 q^-1 + ... + pn q^-n = K*(1 - p10 q-^1) * ... * (1 - pn0 q-^1)
%The parametervector val of polynomial P is defined as. A zero at infinity corresponds to a delay.
%Factorized description: [K p10 ... pn0].
%Standard (defactorized) description: [p0 p1 ... pn].
%
% factorize(obj,poly,ix)
% obj [NsfPolyModel]: The Model object
% poly [opt. char]: Character that defines the polynomial that should be factorized ('A' 'B' ... 'F'). If no argument is
%                   supplied then all polynomials will be factorized
% ix [opt. 2 x 1 pos. int]: Element of poly that should factorized (poly need to be supplied).
%
% EXAMPLE (Factorize SIMO ARMAX Model):
% M = idModels.NsfPolyModel([2 2; 1 2],[2 2]',[2 4]',[2 2]);
% M.A(1,1).val = conv([1 -.9+.1i],[1 -.9-.1i]);
% M.A(1,1).free = [0 0 0];
% M.A(1,2).val(2:end) = [-.6 -.2];
% M.factorize('A');
% M.B(1).val(end) = 2;
% M.factorize('B');

uf = obj.UpdateFlag;
obj.UpdateFlag = false;
if nargin==1 || isempty(polystr)
    obj.factorize({'A' 'B' 'C' 'D' 'E' 'F'});
    return;
end

if iscell(polystr)
    for nP = 1:length(polystr)
        obj.factorize(polystr{nP});
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
        if ~P(r,c).factorized
            if length(P(r,c).val)>21
                warning('The polynomial %s(%i,%i) is of high order (>20)! Results may become inaccurate!',polystr,r,c);
            end
            if length(P(r,c).free)>nk(r,c) && P(r,c).free(nk(r,c)+1); Kfree = true; else; Kfree = false; end
            if ~all(P(r,c).free(nk(r,c)+2:end)) && ~all(~P(r,c).free)
                warning('Unique transformation of %s(%i,%i).free is not possible. All poles will be free after factorization!',polystr,r,c);
            end
            P(r,c).val = idModels.util.factorizePoly(P(r,c).val);
            if ~all(~P(r,c).free) 
                P(r,c).free = [Kfree false(1,nk(r,c)) true(1,length(P(r,c).val)-nk(r,c)-1)];
            end   
            P(r,c).factorized = true;
            P(r,c).std = NaN(1,length(P(r,c).val));
        else
            warning('The Polynomial %s(%i,%i) is factorized already! -> Skipping!',polystr,r,c);
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
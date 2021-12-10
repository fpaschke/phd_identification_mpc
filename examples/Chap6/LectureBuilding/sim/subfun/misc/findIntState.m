function ix = findIntState(M,tol,type)
%FINDINTSTATE Find stateindex of input integration.
%for Statespace object M. M needs to be in modal form!
%If type 
% Function tries to find stateindex by checking the following criteria:
%   1.) diagonal element of A that is 1 and all elements in that row and 
%       column are 0 (up to tolerance tol). This crition should find
%       integrator state indices.
%   2.) all elements of B that correspond to deterministic inputs need to be 
%       0 as well (up to tol).
%
%ix = findDriftState(M,tol)
%   ix [int vector]:        possible state indices.
%   M [matlab Ss object]:   State space object. Needs to be in Modal form!
%   tol [pos. double]:      Tolerance.

if nargin < 2; tol = 1e-8; end
if nargin < 3; type = 'noise'; end
assert(tol>=0,'tol needs to be >=0');
ix = find(abs(diag(M.A)-1)<tol);
isIntState = false(length(ix),1);
for n = 1:length(ix)
    colA = M.A(:,ix(n));
    colA(ix(n)) = NaN;
    colMaxA = max(abs(colA));
    
    rowA = M.A(ix(n),:);
    rowA(ix(n)) = NaN;
    rowMaxA = max(abs(rowA));
    
    if strcmpi(type,'noise') % If pure noise integration the row in B needs to be 0
        if isa(M,'ss') && isfield(M.InputGroup,'Measured')
            rowMaxB = max(M.B(ix(n),M.InputGroup.Measured));
        else
            rowMaxB = max(M.B(ix(n),:)); 
        end
    elseif strcmpi(type,'input') % Only one entry is one
        if sum(abs(M.B(ix(n),:)-1) < tol) == 1
        	rowMaxB = 0;
        else
            rowMaxB = 2*tol;
        end
    end
    
    if colMaxA < tol && rowMaxA < tol && rowMaxB < tol
        isIntState(n) = true;
    else
        isIntState(n) = false;
    end
end
ix = ix(isIntState);
end


function varargout = polyStruct2cell(varargin)
%Converts polynomial structures P1, P2, ... to cell.
%
% [P1c, ..., Pnc] = polyStruct2cell(P1,P2,...,Pn,varargin)
% P1,...,Pn [struct array]: Polynomial structures with fields val, free and factorized.
% P1c,...,Pnc [cell array]: Cell array containing coefficient values of P1 ... Pn
% varargin: Optional name value pairs:
%   SetNanTo0 [logical]: NaN Values will be set to 0 if true (Def. false).
%   DefactorizePoly [logical]: If polynomial struct is in factorized form, then it will be defactorized.
%   ImagTol [pos]: If polynomial gets defactorized then ignore imaginary parts that are less or equal then ImagTol


nP = find(cellfun(@(xi) isstruct(xi),varargin),1,'last');
p = inputParser;
addParameter(p,'SetNanTo0',false,@(x) x==1 || x==0);
addParameter(p,'DefactorizePoly',true,@(x) x==1 || x==0);
addParameter(p,'ImagTol',1e-12,@(x) isscalar(x) && x>=0 && isreal(x));
parse(p,varargin{nP+1:end});

varargout = cell(1,nP);
for np = 1:nP
    [r,c] = size(varargin{np});
    varargout{np} = cell(r,c);
    for nr = 1:r
        for nc = 1:c
            if p.Results.DefactorizePoly && isfield(varargin{np}(nr,nc),'factorized') && varargin{np}(nr,nc).factorized
                varargout{np}{nr,nc} = idModels.util.defactorizePoly(varargin{np}(nr,nc).val,p.Results.ImagTol);
            else
                varargout{np}{nr,nc} = varargin{np}(nr,nc).val;
            end
            if p.Results.SetNanTo0
                varargout{np}{nr,nc}(isnan(varargout{np}{nr,nc})) = 0;
            end
        end
    end
end
end

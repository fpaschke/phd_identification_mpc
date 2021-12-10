function [A, B, C, D, E, var_e, eps] = plr(y,u,na,varargin)
%Estimates Polynomial Model 
%
%   A(q^-1) y[k] = E^-1(q)*B(q^-1) u[k] + C(q^-1)./D(q^-1) e[k]
%
%using PLR Algorithm. E is assumed to be diagonal.
%
% [A, B, C, D, E, var_e, e] = plr(y,u,na,nb,nk,nc,nd,ne,varargin)
% y [Ns x ny double or cell of Ns x ny double]: Output data
% u [Ns x nu double or cell of Ns x nu double]: Input data
% na [ny x ny pos int/polystruct]: Degree of A polynomials or polystruct -> see idModels.util.polyStruct()
% nb [ny x nu pos int/polystruct]: Degree of B polynomials or polystruct -> see idModels.util.polyStruct()
% nk [ny x nu pos int]: Input delays (only available if nb is pos int matrix)
% na [ny x ny pos int/polystruct]: Degree of A polynomials or polystruct -> see idModels.util.polyStruct()
% na [ny x ny pos int/polystruct]: Degree of A polynomials or polystruct -> see idModels.util.polyStruct()
% varargin [Name-Value-Pairs or options struct]: see plrOptions() for available options
% Output:
% A [ny x ny cell of double vectors]: Polynomial coefficients of A
% B [ny x nu cell of double vectors]: Polynomial coefficients of B
% vare [ny x ny double]: Covariance of Noise e 
% e [Ns x ny double or cell of Ns x ny double]: Residuals


if ~iscell(y); y = {y}; end
if ~iscell(u); u = {u}; end
nu = size(u{1},2);
ny = size(y{1},2);

% Get Options Struct
ix_opt = find(cellfun(@ischar,varargin),1,'first'); % Find Options in varargin
if ~isempty(ix_opt)  
    opt = idModels.alg.ls.plrOptions(varargin{ix_opt:end}); 
elseif isstruct(varargin{end})
    opt = varargin{end};
    ix_opt = length(varargin);
else
    opt = idModels.alg.ls.plrOptions(); 
    ix_opt = length(varargin)+1;
end

% Extract na, ..., ne
if isempty(na); na = zeros(ny,ny); end
if ix_opt>1 && ~isempty(varargin{1}); nb = varargin{1}; elseif nu == 0; nb = zeros(ny,nu); else; error('You need to specify nb and nk!'); end
if ix_opt>2 && (~isempty(varargin{2}) || isstruct((varargin{1}))); nk = varargin{2}; elseif nu == 0; nk = zeros(ny,nu); else; error('You need to specify nb and nk!'); end
if ix_opt>3 && ~isempty(varargin{3}); nc = varargin{3}; else; nc = zeros(ny,ny); end
if ix_opt>4 && ~isempty(varargin{4}); nd = varargin{4}; else; nd = zeros(ny,ny); end
if ix_opt>5 && ~isempty(varargin{5}); ne = varargin{5}; else; ne = zeros(ny,ny); end

% Get PolyStructs
if isvector(nc) ncvec = true; else; ncvec = false; end
if isvector(nd) ndvec = true; else; ndvec = false; end
if isvector(ne) nevec = true; else; nevec = false; end
a = idModels.util.polyStruct(na,'A'); na = arrayfun(@(ai) length(ai.val)-1,a);
b = idModels.util.polyStruct(nb,'B',nk); nb = arrayfun(@(bi) length(bi.val)-1,b);
c = idModels.util.polyStruct(nc,'A'); nc = arrayfun(@(ci) length(ci.val)-1,c);    
d = idModels.util.polyStruct(nd,'A'); nd = arrayfun(@(di) length(di.val)-1,d);
e = idModels.util.polyStruct(ne,'A'); ne = arrayfun(@(ei) length(ei.val)-1,e);

assert(isequal(size(a),[ny ny]) && (isequal(size(b),[ny nu])) || isempty(b) && nu == 0 ...
    && isequal(size(c),[ny ny])  && isequal(size(d),[ny ny])  && ...
    isequal(size(e),[ny ny]) || isempty(e) && nu == 0)

Ns = length(y);
Nss = cellfun(@length,y);

%% Prefiltering
opt.Prefilter = idModels.util.validatePrefilter(opt.Prefilter,ny);
dmax_pf = max(cellfun(@(wi) length(wi{1})-1,opt.Prefilter));
dmax_y = max(na(:)); % v = (I-D)v + (C-I)e + e
dmax_v = max([nc(:); nd(:)]); % v = (I-D)v + (C-I)e + e
dmax_w = max([max([ne(:)]) - max([nb(:)]) 0]);
dmax = dmax_pf + max([dmax_v dmax_w dmax_y]);

%% Init Outs
A = cell(opt.Iter,1); B = A; C = A; D = A; E = A; var_e = A;
[A{1}, B{1}, C{1}, D{1}, E{1}] = idModels.util.polyStruct2cell(a,b,c,d,e,'SetNanTo0',true);

%% STEP 1: ARX Model für initialisierung
if any(na(:)>0) || any(ne(:)>0) || any(nb(:)>0) 
    if nu == 0 || strcmpi(opt.Init,'ls') || all(nc(:)==0) && all(nd(:)==0) && all(ne(:)==0) 
        if any(na(:)>0) %AR... Model
            [A{1},B{1},var_e{1},eps] = idModels.alg.ls.arx_ls(y,u,a,b,[],'Prefilter',opt.Prefilter); 
        else %OE 
            [E{1},B{1},var_e{1},eps] = idModels.alg.ls.arx_ls(y,u,e,b,[],'Prefilter',opt.Prefilter); 
        end
    else
        if any(na(:)>0) %AR... Model
            [A{1},B{1},var_e{1},eps] = idModels.alg.ls.oe_iv(y,u,a,b); 
        else %OE 
            [E{1},B{1},var_e{1},eps] = idModels.alg.ls.oe_iv(y,u,e,b); 
        end
    end
else
    eps = y;
end

% ARX CASE
if all(nc(:)==0) && all(nd(:)==0) && all(ne(:)==0)
    A = A{1}; B = B{1}; C = C{1}; D = D{1}; E = E{1}; var_e = var_e(1); e = E;
    return;
end

%Prepare Polynomials
for no1 = 1:ny
    for no2 = 1:ny 
        % set leading coeffs to 0 (because code idModels.alg.ls.arx_ls(y,uevw,a,[b c d e]); 
        if ~isempty(c) c(no1,no2).val(1) = 0; end
     	if ~isempty(d) d(no1,no2).val(1) = 0; end
        if ~isempty(e) e(no1,no2).val(1) = 0; end
    end
end 

%% STEP 2: Iter
opt_arxls = idModels.alg.ls.arx_lsOptions('Prefilter',opt.Prefilter);
opt_arxls.PlrFlag = true;
uevw = cell(Ns,1); 
for k = 2:opt.Iter
    % Ay = 1/E(Bu) + C/De = w + v <-> y = (1-A)y + w + v 
    % PLR - Darstellung: Ay = Bu + (C-1)e + (1-D)v + (1-E)w
    for ns = 1:Ns
        w = NaN(Nss(ns),ny); v = w; ee = w;
        for no = 1:ny
            % 1.) Comp. w
            if nu > 0
                w(:,no) = sum(cell2mat(arrayfun(@(ni) filter(B{k-1}{no,ni},E{k-1}{no,no},u{ns}(:,ni)),1:nu,'UniformOutput',false)),2);
            else 
                w(:,no) = zeros(Nss(ns),1);
            end
            % 2.) Comp. v = Ay - w
            if any(nd(:)) > 0
                v(:,no) = sum(cell2mat(arrayfun(@(no2) filter(A{k-1}{no,no2},1,y{ns}(:,no2)),1:ny,'UniformOutput',false)),2); %Ay
            	if nu > 0 % Do not subtract zeros
                    v(:,no) = v(:,no) - w(:,no);
                end
            else
                v(:,no) = zeros(Nss(ns),1);
            end
            % 3.) Comp e
            ee(:,no) = eps{ns}(:,no);
            ee(isnan(ee(:,no)),no) = 0;
%             if ~isequal(opt.Prefilter{no},{1 1}) % Prefilter if necessary
%             	ee(:,no) = filter(opt.Prefilter{no}{2},opt.Prefilter{no}{1},ee(:,no));
%             end
        end
        uevw{ns} = horzcat(u{ns},ee,-v,-w);
    end    
    [A{k},BB,var_e{k},eps] = idModels.alg.ls.arx_ls(y,uevw,a,[b c d e],[],opt_arxls);
    
    % extract Parameters
    B{k} = BB(:,1:nu);
    C{k} = BB(:,nu+1:nu+ny);
	D{k} = BB(:,nu+ny+1:nu+2*ny);
  	E{k} = BB(:,nu+2*ny+1:end);
   	for no = 1:ny
        C{k}{no,no}(1) = 1;
        D{k}{no,no}(1) = 1;
     	E{k}{no,no}(1) = 1;
    end 
end

%% Plot Coefficients
% Plotte Koeffizienten der Polynome A B C
if opt.Plot
    idModels.alg.ls.plotCoeff(A,B,C,D,E);
end

%% Output
A = A{end}; 
B = B{end}; 
C = C{end}; 
if ncvec; C = arrayfun(@(k) C{k,k},(1:ny)','UniformOutput',false);  end
D = D{end};
if ndvec; D = arrayfun(@(k) D{k,k},(1:ny)','UniformOutput',false);  end
E = E{end};
if nevec; E = arrayfun(@(k) E{k,k},(1:ny)','UniformOutput',false);  end
var_e = var_e{end};
return;
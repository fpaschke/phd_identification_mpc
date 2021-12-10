function [A,B,vare,e] = arx_ls(y,u,na,nb,nd,varargin)
%Estimate MIMO ARX Model A y(k) = q^-nd B u(k) + e(k) with ny outputs and nu inputs
%with Least squares Method.
%
% [A,B,vare,e] = arx_ls(y,u,na,nb,nd,varargin)
% y [Ns x ny double or cell of Ns x ny double]: Output data
% u [Ns x nu double or cell of Ns x nu double]: Input data
% na [ny x ny pos int/polystruct]: Degree of A polynomials or polystruct -> see idModels.util.polyStruct()
% nb [ny x nu pos int/polystruct]: Degree of B polynomials or polystruct -> see idModels.util.polyStruct()
% nd [ny x nu pos int]: Input delays (only available if nb is pos int matrix)
% varargin [Name-Value-Pairs or options struct]: see arx_lsOptions() for available options
% Output:
% A [ny x ny cell of double vectors]: Polynomial coefficients of A
% B [ny x nu cell of double vectors]: Polynomial coefficients of B
% vare [ny x ny double]: Covariance of Noise e 
% e [Ns x ny double or cell of Ns x ny double]: Residuals
%
% EXAMPLES: see +idModels\+alg\+ls\test folder

% in Cell Konvertieren
flag = 0;
if ~iscell(y) 
    y = {y};
    flag = 1;
end
if ~iscell(u)
    u = {u};
end

if length(varargin) > 1 || isempty(varargin) 
    opt = idModels.alg.ls.arx_lsOptions(varargin{:}); 
elseif isstruct(varargin{1}) && length(varargin) == 1
    opt = varargin{1};
else
   error('Invalid Input Argument!') 
end
if isfield(opt,'PlrFlag'); plrflag = opt.PlrFlag; else; plrflag = false; end

if ~isstruct(na) 
    a = idModels.util.polyStruct(na,'A');
else
    a = na;
    na = arrayfun(@(ai) length(ai.val)-1,a);
end
if ~isstruct(nb)
    b = idModels.util.polyStruct(nb,'B',nd);
else
    b = nb;
    nb = arrayfun(@(bi) length(bi.val)-1,b);
end

nu = length(nb);
ny = length(na);
if nu == 0
   B = cell(ny,0); 
end

use_iv = 0;
if ~isempty(opt.Instruments)
    z = opt.Instruments;
    use_iv = 1; % Hilfsvariablen nutzen
    if ~iscell(z)
        z = {z};
    end
end

opt.Prefilter = idModels.util.validatePrefilter(opt.Prefilter,ny);
%assert(ny== 1 || all(cellfun(@(pfi) isequal(pfi,{1 1}),opt.Prefilter)),'Prefilters arent supported for SIMO/MIMO Systems yet!');
dmax_w = cellfun(@(wi) length(wi{1})-1,opt.Prefilter);
dmax = max(max([na nb])) + max(dmax_w);

Y = []; PHI = []; ZZ = [];
for k = 1:length(y)
    N = size(y{k},1);
    nr = N - dmax;
    Yk = cell(ny,1);
    for no = 1:ny
        yp = []; up = []; zp = [];
        % Preproc Input 
        Ap = opt.Prefilter{no}{2}; Bp = opt.Prefilter{no}{1};
        if ~isequal(Bp,1) || ~isequal(Ap,1)
            nmax = max([length(Bp) length(Ap)]-1);
            % NOTE: Die Berechnung der ICs ist sinnvoll für HP Vorfilterung 
            % für TP Vorfilterung sollte zeros(1,nmax) durch y{k}(1,no2)*ones(1,nmax) ersetzt werden
            for no2 = 1:ny
                yp(:,no2) = filter(Bp,Ap,y{k}(:,no2),filtic(Bp,Ap,zeros(1,nmax),y{k}(1,no2)*ones(1,nmax)));
             	if use_iv
                    %zp(:,no2) = filter(Bp,Ap,z{k}(:,no2),filtic(Bp,Ap,zeros(1,nmax),z{k}(1,no2)*ones(1,nmax)));
                    zp(:,no2) = z{k}(:,no2); % prefiltering of instruments is handled outside of arx_ls
                end
            end
            for ni = 1:nu 
                if ni <= nu - plrflag*3*ny
                    up(:,ni) = filter(Bp,Ap,u{k}(:,ni),filtic(Bp,Ap,zeros(1,nmax),u{k}(1,ni)*ones(1,nmax)));
                else %The last columns are the sequences e, v ,w. They shouldnt be filterd. plr function.
                    up(:,ni) = u{k}(:,ni);
                end
            end
        else
            yp = y{k};
            up = u{k};
            if use_iv
                zp = z{k};
            end
        end
        
    	% Y
        Yk{no} = yp(dmax+1:N,no);
        
        % A
        Phi_yk = []; Z_yk = [];
        for no2 = 1:ny
            ix_c = a(no,no2).val ~= 0 | a(no,no2).free; % only care about parameters which are not 0 or free
            Phi = zeros(nr,sum(ix_c)-1); %Phi = zeros(nr,na(no,no2));
            if use_iv
                Z = zeros(nr,sum(ix_c)-1); %Z = zeros(nr,na(no,no2));
            end
            nc = 1;
            for i=find(ix_c(2:end)) %1:na(no,no2)
                Phi(:,nc) = - yp(dmax+1-i:N-i,no2); nc = nc + 1;
                if use_iv
                    Z(:,nc) = - zp(dmax+1-i:N-i,no2);
                end
            end
            % Subrahiere festen Anteil
            ix_s = [false ~a(no,no2).free(2:end) & a(no,no2).val(2:end)~=0]; % indizes der festen nicht 0 parameter 
            if any(~a(no,no2).free([false ix_c(2:end)]))
                Yk{no} = Yk{no} - Phi(:,~a(no,no2).free([false ix_c(2:end)]))*a(no,no2).val(ix_s)';
                %Yk{no} = Yk{no} - Phi(:,~a(no,no2).free(2:end))*a(no,no2).val([false ~a(no,no2).free(2:end)])';
            end
            Phi_yk = [Phi_yk Phi(:,a(no,no2).free([false ix_c(2:end)]))]; %Phi_yk = [Phi_yk Phi(:,a(no,no2).free(2:end))];
            if use_iv
                Z_yk = [Z_yk Z]; %Z_yk = [Z_yk Z(:,a(no,no2).free(2:end))];
            end
        end

        % Die mit den Eingangswerten
        Phi_uk = [];
        for ni = 1:nu
            ix_c = b(no,ni).val ~= 0 | b(no,ni).free; % only care about parameters which are not 0 or free
            Phi = zeros(nr,sum(ix_c)); %Phi = zeros(nr,nb(no,ni)+1);
            nc = 1; % Column index
            for j = find(ix_c)-1 %0:nb(no,ni)
                Phi(:,nc) = up(dmax+1-j:N-j,ni); nc = nc + 1;
            end
            % Subrahiere festen Anteil
            ix_s = ~b(no,ni).free & b(no,ni).val~=0; % indizes der festen nicht 0 parameter 
            if any(~b(no,ni).free(ix_c))
                Yk{no} = Yk{no} - Phi(:,~b(no,ni).free(ix_c))*b(no,ni).val(ix_s)';
                %Yk{no} = Yk{no} - Phi(:,~b(no,ni).free)*b(no,ni).val(~b(no,ni).free)';
            end
            Phi_uk = [Phi_uk Phi(:,b(no,ni).free(ix_c))];
        end
        Phi_k{no} = [Phi_yk Phi_uk]; 
        if use_iv
        	Z_k{no} = [Z_yk Phi_uk]; 
        end
    end
    
 	% cct
    Y = [Y; cell2mat(Yk)];
    PHI = [PHI; 
           blkdiag(Phi_k{:})]; 
	if use_iv
        ZZ = [  ZZ; 
                blkdiag(Z_k{:})]; 
    end
end

% Die MKQ Lösung und Residuen berechnen
if ~use_iv
	P = PHI\Y; % [a11,1 ... a11,na a12,1 ... a12,na]
else
 	P = (ZZ'*PHI)\(ZZ'*Y);
end

% Parameter extrahieren
ixs = 1;
for no1 = 1:ny
    for no2 = 1:ny
        ixe = ixs + sum(a(no1,no2).free) - 1;
        A{no1,no2} = a(no1,no2).val;
        A{no1,no2}(:,a(no1,no2).free) = P(ixs:ixe)';
        ixs = ixe + 1;
    end

    for ni = 1:nu
        ixe = ixs + sum(b(no1,ni).free) - 1;
        B{no1,ni} = b(no1,ni).val;
        B{no1,ni}(:,b(no1,ni).free) = P(ixs:ixe)';
        ixs = ixe + 1;
    end
end

% Residuen zu erzeugen
if nargout > 2
	E = Y - PHI*P; 
  	ns = 1;
    e = arrayfun(@(ns) zeros(0,ny),1:length(y),'UniformOutput',0)';
    for j = 1:length(y)
        if size(y{j},1)-dmax>0
            ne = ns + ny*(size(y{j},1)-dmax) - 1;
            e{j,1} = [  NaN(dmax,ny); 
                        reshape(E(ns:ne),size(y{j},1)-dmax,ny)];
            ns = ne + 1;
        end
    end
    E = cell2mat(e);
    E = E(~isnan(E(:,1)),:);
    if size(E,1)-length(P) > 0
        vare = E'*E/(length(E)-length(P));
    else
        vare = E'*E/length(E);
    end
 	if flag 
        e = e{1};
    end
end
return
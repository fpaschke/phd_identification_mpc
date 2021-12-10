function [A,B,var_e,e] = oe_iv(y,u,na,nb,nd,varargin)
%Estimate MIMO OE Modell
%  y(k) = B./F u(k) +  v(k)
%using Instrumental variable method.
%
% [A,E,varv,v] = oe_iv(y,u,nf,nb,nd,varargin)
% y [Ns x ny double or cell of Ns x ny double]: Output data
% u [Ns x nu double or cell of Ns x nu double]: Input data
% nf [ny x nu pos int/polystruct]: Degree of F polynomials or polystruct -> see idModels.util.polyStruct()
% nb [ny x nu pos int/polystruct]: Degree of B polynomials or polystruct -> see idModels.util.polyStruct()
% nd [ny x nu pos int]: Input delays (only available if nb is pos int matrix)
% varargin [Name-Value-Pairs or options struct]: see oe_ivOptions() for available options
% Output:
% F [ny x nu cell of double vectors]: Polynomial coefficients of F
% B [ny x nu cell of double vectors]: Polynomial coefficients of B
% varv [ny x ny double]: Covariance of Noise v
% v [Ns x ny double or cell of Ns x ny double]: Residuals
%
% EXAMPLES: see +idModels\+alg\+ls\test folder

if length(varargin) >= 2 || isempty(varargin) 
    opt = idModels.alg.ls.oe_ivOptions(varargin{:}); 
elseif isstruct(varargin{1}) && length(varargin) == 1
    opt = varargin{1};
else
   error('Invalid Input Argument!') 
end

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

% in Cell Konvertieren
flag = 0;
if ~iscell(y) 
    y = {y};
    flag = 1;
end
if ~iscell(u)
    u = {u};
end

%Init
Ne = length(y);
Ns = cellfun(@length,y);
Nmax = max(Ns);
if isequal(opt.Ls_init_samples,Inf)
    Ninit = max(Ns);
elseif strcmpi(opt.Ls_init_samples,'auto')
    Ninit = 5*na;
else
    Ninit = opt.Ls_init_samples;
end
opt.Prefilter = idModels.util.validatePrefilter(opt.Prefilter,ny);

z = []; % first estimate simple arx model with ls technique
A = cell(opt.Iter,1); B = A; var_e = A;
for j = 1:opt.Iter
    [A{j},B{j},var_e{j},e] = idModels.alg.ls.arx_ls(y,u,a,b,[],z,'Prefilter',opt.Prefilter);
	
    %% Filter: Simulate System
    % get sequences W*A*y = W*B*u = w
    if ny == 1
        Ai_num = {1}; Ai_den = A{j};
        if opt.EstIc>0
            PSI = idModels.alg.lti.getPsi(A{j}{1},Nmax); %we need the full PSI because we use it for simulation e = e + PSI*e0 
        end
    else
        Ga = tf(idpoly(A{j},[],num2cell(ones(ny,1)),[],[],1),'Noise'); % Compute A^-1
        Ai_num = Ga.num; Ai_den = Ga.den;
        if opt.EstIc>0
            PSI = idModels.alg.lti.getPsi(Ga,Nmax); %we need the full PSI because we use it for simulation e = e + PSI*e0
        end
    end
    %P = idModels.alg.lti.stabilizePoly(Ai_den{1,1}); %detA;
    %for no1 = 1:ny for no2 = 1:ny Ai_den{no1,no2} = P; end; end;
        
  	for ns = 1:Ne
        % w = W*B*u
        w = [];
        for no = 1:ny 
            w(:,no) = sum(cell2mat(arrayfun(@(ix_u) filter(conv(opt.Prefilter{no}{1},B{j}{no,ix_u}),opt.Prefilter{no}{2},u{ns}(:,ix_u)),1:nu,'UniformOutput',0)),2); % get w = sum(Bi*ui)
        end
      	% calc forced output y = A^(-1)*w
        for no = 1:ny
            z{ns}(:,no) = sum(cell2mat(arrayfun(@(ix_w) filter(Ai_num{no,ix_w},Ai_den{no,ix_w},w(:,ix_w)),1:ny,'UniformOutput',0)),2); % get w = sum(Bi*ui)
        end
        if opt.EstIc>0
            Ni = min(Ninit,Ns(ns));
            yp = cell2mat(arrayfun(@(no) filter(opt.Prefilter{no}{1},opt.Prefilter{no}{2},y{ns}(1:Ni,no)),1:ny,'UniformOutput',0));
            D = (yp-z{ns}(1:Ni,:))';
            y0 = PSI(1:Ni*ny,:)\D(:); % calc init cond
            z0 = PSI(1:ny*Ns(ns),:)*y0;
            z{ns} = z{ns} + reshape(z0,ny,Ns(ns))';
        end
        %close all; plot([z{ns} y{ns}]);
    end
end

%% Plot Coefficients
% Plotte Koeffizienten der Polynome A B C
if opt.Plot
    idModels.alg.ls.plotCoeff([],B,[],[],A);
end

%% Output
A=A{end};B = B{end}; var_e = var_e{end};
if flag 
    e = e{1};
end
function [parTable, misc] = sweepPar(varyingpars,estfun,varargin)
% Sweeps parameters eg for model identification and returns result table.
%
% [parTable, misc] = sweepPar(varyingpars,estfun,varargin)
% parTable [table]: Resulttable. First column contains costs. The other columns contain the different Parametercombinations
% varyingpars [string/cell of strings]: Parameters in varargin to be varied
% estfun [function handle]: Function handle that computes costs.
%                           Output arguments are costs meaning that they need to be a scalar double. The parTable will be sorted by
%                           by the first cost. estFun needs to be able to handle the arguments in varargin, i. e. it will be
%                           called by estfun(varargin{:})
% varargin: The Arguments to be passed to estfun in the specified order including the varyingpars, i e varargin needs to
%           contain varyingpars as name value pair arguments (See exaple below).
%
% EXAMPLE:
% rng(0); u = kron(randn(20,1),ones(50,1)); e = .1*randn(length(u),1);
% A = conv([1 -.9],[1 -.8]); B = [0 1]; C = [1 -.9];
% y = filter(B,A,u) + filter(C,A,e);
% res = idModels.test.sweepPar({'Order' 'Type'},@idModels.NsfPolyModel.estModel,y,u,'Order',1:5,'Type',{'ARX' 'ARMAX' 'BJ'},'HpVal',10);
% [mat, lab] = tab2mat(res,1); bar(sqrt(mat)); ylim(sqrt([2e-2 8e-2])); xlabel('Order of Polynomials');
% ylabel('RMMSE of 10-Step-Prediction-Errors'); legend(lab{2});

ix_costfun = find(cellfun(@(x) isa(x,'function_handle'),varargin),1,'first');
if ischar(varyingpars)
    varyingpars = {varyingpars};
end

p = inputParser;
p.KeepUnmatched = true;
addRequired(p,'estfun',@(x) isa(x,'function_handle'));
addRequired(p,'varyingpars',@iscell);
parse(p,estfun,varyingpars)

%find varying parameters
parn = p.Results.varyingpars;
ix_vp = cellfun(@(vp) find(cellfun(@(x) isequal(x,vp),varargin)),varyingpars) + 1;
vals_vp = varargin(ix_vp);

% get max number of variations for each parameter
nPar = length(parn);
ix_max = []; %maximum values 
iscellstrpar = []; nVecPar = [];

for np = 1:nPar
    if isnumeric(vals_vp{np}) || (iscell(vals_vp{np}) && ischar(vals_vp{np}{1}))
        ix_max = [ix_max length(vals_vp{np})]; % maximum possible values of parameters
        if isvector(vals_vp{np}) % Make sure to be row vector
            vals_vp{np} = vals_vp{np}(:);
        end
        if iscell(vals_vp{np}) && ischar(vals_vp{np}{1})
            iscellstrpar = [iscellstrpar 1];
        else
            iscellstrpar = [iscellstrpar 0];
        end
        nVecPar = [nVecPar 1];
    elseif  iscell(vals_vp{np}) && isnumeric(vals_vp{np}{1})
        iscellstrpar = [iscellstrpar 0];
        nVecPar(np) = length(vals_vp{np});
        for nvp = 1:nVecPar(np)
            ix_max = [ix_max length(vals_vp{np}{nvp})];
        end
    else
        error('Check format of varying parameters!');
    end    
end
nmods = prod(ix_max); 

% cost 
Nout = nargout(estfun);
efc_out = cell(Nout,1);
%cost = NaN(nmods,Ncost);

% Init Parametertable
parTable = table();

% loop through all possible combinations
ix_act = ones(size(ix_max)); %saves state of actual iteration
n_iter = 1;
h = waitbar(0,'Calculating...');
while 1
    % build varargin argument to pass to estfun
	str = []; P = cell(1,nPar);
    for np = 1:nPar 
        if nVecPar(np)>1 % case: vector Parameter
            val = [];
            for nvp = 1:nVecPar(np)
            	ii = sum(nVecPar(1:np-1)) + nvp;
                val =  [val vals_vp{np}{nvp}(ix_act(ii))];
            end
            varargin{ix_vp(np)} = val;
            actval = mat2str(val);
            P{np} = val;
        else %case: string par or scalar par
            ii = sum(nVecPar(1:np));
            if ~iscellstrpar(np)
                varargin{ix_vp(np)} = vals_vp{np}(ix_act(ii),:);
                actval = num2str(varargin{ix_vp(np)});
            else
                varargin{ix_vp(np)} = vals_vp{np}{ix_act(ii)};
                actval = varargin{ix_vp(np)};
            end
            P{np} = vals_vp{np}(ix_act(ii),:);
        end
        if np < nPar
            str = [str varargin{ix_vp(np)-1} '=' actval '; '];
        else
            str = [str varargin{ix_vp(np)-1} '=' actval];
        end
    end
    
	
	% Check for Cancel button press
    if getappdata(h,'canceling')
        delete(h);
        return;
    end
    % output state
    waitbar(n_iter/nmods,h,['Parameters: ' str]);
    fprintf(['Parameters: ' str '\n']);
    
	% now call estfun and save outputs
    ix_tmp = mat2cell(ix_act,1,ones(1,length(ix_act)));
    try
        [efc_out{:}] = estfun(varargin{:});
        cost(n_iter,:) = efc_out{1};
        if Nout > 1
            misc(n_iter,1) = efc_out(2:end);
        end
    catch ME
        warning(['Some Error Occured during estimation! -> Skipping Modelstructure: ' str]);
        disp(ME.message)
        mdls{n_iter,1} = [];
        cost(n_iter,1) = Inf;
    end
    
    % write parametercombination to table
    parTable = [parTable; {cost(n_iter,:)} P];
    
    n_iter = n_iter + 1; 
    % abort if all ix_act have reached thier maximum
    if all(ix_act == ix_max) 
        break;
    end
    
	% increment is_act
    ix = find(ix_act<ix_max,1,'first'); % find first parameter which didnt reach his limit
    ix_act(ix) = ix_act(ix) + 1; % now increment ix_act
    ix_act(1:(ix-1)) = 1; % now reset all previous parameters to 1  
end
parTable.Properties.VariableNames = ['Cost' varyingpars];
pTab = parTable; % Sicherstellen, dass nach dem ersten Element sortiert wird
pTab.Cost = pTab.Cost(:,1);
[~,ix] = sortrows(pTab,'Cost');
parTable = parTable(ix,:);
if Nout > 1
    misc = misc(ix);
end
delete(h);
end
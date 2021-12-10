function [mat,vals] = tab2mat(tab,ixcost)
%TAB2MAT Sort Table to wrt to ixcost and output an matrix.

if nargin == 1
    ixcost = 1;
end

Cost = tab.(tab.Properties.VariableNames{1});
varnames = tab.Properties.VariableNames(2:end);
nvar = length(varnames );

nr = size(Cost,1);
tabvar = cell(nr,0);
for nv = 1:nvar
    if ~iscell(tab.(varnames{nv}))
        vals{nv} = num2cell(sort(unique(tab.(varnames{nv}))));
        col = num2cell(tab.(varnames{nv}));
    else
        vals{nv} = sort(unique(tab.(varnames{nv})));
        col = tab.(varnames{nv});
    end
    tabvar = [tabvar col];
end
ix = ones(1,nv);
ixe = cellfun(@length,vals);

if length(ixe)==1
    mat = NaN(ixe,1);
else
    mat = NaN(ixe);
end
while 1
    var_val = arrayfun(@(nv) vals{nv}{ix(nv)},1:nvar,'UniformOutput',0);
    nrow = find(arrayfun(@(r) isequal(tabvar(r,:),var_val),1:nr));
	idx = num2cell(ix);
    mat(idx{:}) = Cost(nrow,ixcost);
    if all(ix == ixe)
       break; 
    end
  	idx = find(ix<ixe,1,'first'); % find first parameter which didnt reach his limit
    ix(idx) = ix(idx) + 1; % now increment ix_act
    ix(1:(idx-1)) = 1; % now reset all previous parameters to 1  
end
end


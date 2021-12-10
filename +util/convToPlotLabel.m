function lab = convToPlotLabel(lab,interpreter)
% Converts a string to a valid nice plotlabel by replacing characters such as '_' by '\_' such that they 
% will displayed nicely. 
%
% clab = convToPlotLabel(lab)
% lab [char or cellstr]: the labels to be converted.
% clab [char or cellstr]: the converted labels.

if nargin == 1
   interpreter = 'latex'; 
end

cellin = 1;
if ~iscell(lab)
    lab = {lab};
    cellin = 0;
end

for k = 1: length(lab)
    if strcmpi(interpreter,'latex')
        lab{k} = strrep(lab{k},'\_','_'); % damit nicht \\_ rauskommt falls \_ schon vorhanden 
        lab{k} = strrep(lab{k},'_','\_');
     	lab{k} = strrep(lab{k},'ä','\"a');
        lab{k} = strrep(lab{k},'ö','\"o');
        lab{k} = strrep(lab{k},'ü','\"u');
        lab{k} = strrep(lab{k},'ß','{\ss}');
        lab{k} = strrep(lab{k},'degC','$^{\circ}$C');
        lab{k} = strrep(lab{k},'deg','$^{\circ}$');
        lab{k} = strrep(lab{k},'%','\%');
        % TO BE CONTINUED HERE
    else
        error('Not implemented yet!');
    end
end

if ~cellin
   lab = lab{1}; 
end
end
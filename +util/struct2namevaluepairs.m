function nvp = struct2namevaluepairs(str)
%STRUCT2NAMEVAL Function converts a structure to name valuepair cell.
% 
% nvp = struct2namevaluepairs(str)
% str [struct]: Structure to be converted to name value pairs.
% nvp [cell]: Cell of name value pairs.
%
% EXAMPLE:
% str.alpha = 90;
% str.beta = 180;
% nvp = struct2namevaluepairs(str)

assert(isstruct(str),'Input needs to be a structure!');
nvp = {};
fnames = fieldnames(str);
for nf = 1:length(fnames)
    nvp = {nvp{:} fnames{nf} str.(fnames{nf})};
end
end


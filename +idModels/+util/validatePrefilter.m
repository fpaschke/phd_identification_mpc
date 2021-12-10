function W = validatePrefilter(W,ny)
%VALIDDATEPREFILTER Summary of this function goes here
%   Detailed explanation goes here
if isempty(W)
    for no = 1:ny  
        W{no} = {1 1};
    end
end
if ny == 1 && length(W) == 2
    W = {W};
end
for no = 1:ny  
    if isempty(W{no}) || isscalar(W{no})
        W{no} = {1 1};
    end
end
end


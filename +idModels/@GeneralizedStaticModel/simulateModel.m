function y = simulateModel(obj,u,opt)
% SIMULATE Simulates  Model for given Inputs u.

import idModels.func.*;

y = NaN(size(u,1),obj.OutputDimension);
for ny = 1:obj.OutputDimension
    fun = str2func(obj.Fun{ny});
    try
        y(:,ny) = fun(u,obj.Parameters{ny});
    catch
        [y(:,ny),~] = fun(u,obj.Parameters{ny}); % in case of deal and fun returning gradient
    end
    y(:,ny) = y(:,ny) + obj.OutputOffset(ny);
end
end
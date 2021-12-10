classdef (Abstract) StaticModel < idModels.IdModel
    %STATICMODEL Defines a abstract class representing static Models of kind y = f(u,p) + e where e is assumed 
    %to be white noise with fixed covariancematrix.
    
    properties
    end
    
    methods %public methods
        function obj = StaticModel(varargin) % Init from  Name Value Pairs
            obj = obj@idModels.IdModel(varargin{:});
        end
    end %public methods
end
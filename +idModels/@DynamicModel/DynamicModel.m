classdef (Abstract) DynamicModel < idModels.IdModel
    %DYNAMICMODEL Represents any kind of discrete time dynamic model.
    
    properties
    end
    
    methods
        function obj = DynamicModel(varargin) % Init from  Name Value Pairs
            obj = obj@idModels.IdModel(varargin{:});
        end
    end
end


classdef SoftplusModel < idModels.GeneralizedStaticModel
    %SOFTPLUSMODEL Defines a sum of multiple Softplus models:
    %             K 
    %           -----
    % y(u,p) =   >   p_k(1)*ln(1 + exp(p_k(2) * ( u - p_k(3) )))   
    %           -----
    %            k=1
    % where p_k is the 3 x 1 parametervector of summand k: 
    % p_k(1)*p_k(2) defines the slope for u->inf
    % p_k(3) is the shifting parameter
    % p_k(2) defines the curvature at p_k(3).
 	% The Number of Summands/Pieces is determined by length of the Parametervector, i.e by length(obj.Parameters{ny}),
    % where ny is the number of the output.
    %
    % NOTE: Only SISO SoftplusModel are defined currently!
    %
    % obj = SoftplusModel(Nin,Nout,varargin)
    % obj [idModels.SoftplusModel]: The object to be created
    % Nin [pos. int]: Number of inputs
    % Nout [pos. int]: Number of outputs
    % varargin [Name-Valuepairs]:   Any valid Property of BModel/StaticModel/GeneralizedStaticModel can be set by Name value
    %                               pairs
    %   
    % EXAMPLES: 
    %
    % Normal Softplus
    % SP1 = idModels.SoftplusModel('Name','SoftPlus','Parameters',[2 1 10],'Free',[1 1 1]);
    % u = (-100:100)';
    % y = 3*(u-20).*(u>=20) + 100 + 5*randn(length(u),1);
    % SP1.identify(y,u,'EstimateOutputOffset',true,'CheckGradient',true);
    % SP1.showModel('-r','Data',{y u})
    % 
    % Twosided Softplus
    % SP2 = idModels.SoftplusModel('Name','SoftPlus','Parameters',[1 -1 5 1 1 30],'Free',[1 1 1 1 1 1],'OutputOffset',10)
    % u = (-20:.1:40)';
    % y = -3*(u-10).*(u<=10) + 3*(u-25).*(u>=25) + 5*randn(length(u),1) + 30;
    % SP2.identify(y,u,'EstimateOutputOffset',true,'CheckGradient',true,'Solver','matlab_lsqnonlin');
    % SP2.showModel('-r','Data',{y u})
    % SP2.show('Autocovariance',y,u,'Significance',.01)
    
    properties (Dependent = true)
        Npieces %[ny x 1 pos. integer]: indicates of how many Pieces/Summands Softplus function consists. 
    end
    
    methods
        %% CONSTRUCTOR
        function obj = SoftplusModel(varargin) % Init from  Name Value Pairs
            obj = obj@idModels.GeneralizedStaticModel(varargin{:});
            assert(obj.InputDimension==1,'Only Models with one input are currently defined');
         	obj.Fun = arrayfun(@(i) 'idModels.func.fun_softplus',(1:obj.OutputDimension)','UniformOutput',false);
        end
        
        %% Getter and setter
        function val = get.Npieces(obj)
            val = cellfun(@(pi) length(pi)/3,obj.Parameters);
        end
    end
    
    %% PROTECTED METHS
    methods (Access = protected)
        function ok = checkParameters(obj)
        % Checks weather Parameters has an valid value
            for no = 1:obj.OutputDimension
                ok(no) = mod(length(obj.Parameters{no}),3) == 0;
            end
            ok = all(ok);
        end
    end
    
 	%% Static Methdos
    methods (Static = true)
        function obj = loadobj(str)
            obj = idModels.SoftplusModel(length(str.InputName),length(str.OutputName));
            obj = obj.set(str);
        end
    end
end


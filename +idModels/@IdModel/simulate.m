function [y, varargout] = simulate(obj,u,varargin)
%SIMULATE Performs simulation of a model for given inputs u.
%This function calls the appropriate simulateModel method for Model simulation.
%See simulateModel function for further information.
%
%[y, varargout] = simulate(obj,u,varargin)
%
% obj [IdModel]: A model object.
% u [(cell of) N x nu double]: Input Data.
% y [(cell of) N x ny double]: Simulated Output Data.
% varargin: specific inputs to be passed to simulateModel (e. g. initial state for dynmaic models). See simulateModel function for further information. 
% varargout: specific outputs returned from simulateModel. See simulateModel function for further information.

%% Extract/Validate Data
if ~iscell(u); cellin = false; else; cellin = true; end
[~,u] = obj.getRawData([],u);
%obj.checkData([],u);

%% Simulate
Ns = length(u);
y = cell(Ns,1);
out = cell(1,nargout-1);
varargout = cell(1,nargout-1);
for ns = 1:Ns
	[y{ns}, out{:}] = obj.simulateModel(u{ns},varargin{:});
    if nargout > 1
        for no = 1:nargout-1
            if ns == 1
                varargout{no} = out(no);
            else
            	varargout{no} = [varargout(no); out(no)];
            end
        end
    end
end
if ~cellin 
    y = y{1};
    varargout = cellfun(@(out) out{1},varargout,'UniformOutput',false);
end
end
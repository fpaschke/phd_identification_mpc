function opt = identifyModelOptions(obj,varargin)
%IDENTIFYOPTIONS Generates options struct for identify function of PolynomialModel.
%
% opt = identifyOptions(varargin)
% obj [PolynomialModel]: The PolynomialModel Object.
% varargin: Name Value Pairs
%   'FirstDerivative' [-1/0/1]:	1/-1 enforces 1st derivative of each inputchannel to be pos/neg (monotonically increasing/decreasing function) 
%                               0 means no restriction (default: 0). The option uses Nconstraints distributed by linspace(umin,umax,Nconstraints)
%                               for each inputchannel to ensure the restrictions on monotonicity
%   'SecondDerivative' [-1/0/1]: 1/-1 enforces 2nd derivative of each inputchannel to be pos/neg (monotonically increasing/decreasing 1st derivative) 
%                               0 means no restriction (default: 0). The option uses Nconstraints distributed by linspace(umin,umax,Nconstraints)
%                               for each inputchannel to ensure the restrictions on 2nd Derivative
%   'Nconstraints' [pos. int]: Number of constraints for each input regarding 'FirstDerivative' and 'SecondDerivative'

assert(obj.OutputDimension == 1,'Identification for MultiOutput Systems is unsupported yet!');
p = inputParser();
addParameter(p,'Nconstraints',10,@(x) mod(x,1)==0 && x >= 0);
addParameter(p,'FirstDerivative',0,@(x) x == 1 || x == -1 || x == 0);
addParameter(p,'SecondDerivative',0,@(x) x == 1 || x == -1 || x == 0);
parse(p,varargin{:});
opt = p.Results;
end
function [Ko,poles,Ak] = kalmanGain(S,varargin)
%KALMANFILTER Returns gain of stationary kalmanfilter from Statespace 
%model in innovations form. The statespace object needs to have noise
%channel, i.e. S.InputGroup.Noise shouldnt be empty.
%
%Ko = kalmanGain(S,varargin)
% Ko [n x ny double]: Kalman gain.
% S [ss]: Matlab Statespace object in innovations form. S.InputGroup.Noise shouldnt be empty.
% varargin: Optional name value pairs
%   'StabilityCheck':   performs stability check of observer and ouputs
%                       warning if observer will be not stable.

p = inputParser();
addParameter(p,'StabilityCheck',true,@(x) x==1 || x==0);
parse(p,varargin{:})

ix_e = S.InputGroup.Noise;
ix_n = any(S.B(:,ix_e),1);
assert(~isempty(ix_e),'The StateSpace object isnt in innovations form!');
Ko = S.B(:,ix_e);
Ko(:,ix_n) = Ko(:,ix_n)./diag(S.D(:,ix_e(ix_n)))';
if p.Results.StabilityCheck || nargout > 1
    Ak = S.A-Ko*S.C;
    poles = eig(full(Ak));
    abs_rts = abs(poles);
    if ~all(abs_rts<1)
        warning('Observer/Predictor is unstable!');
    end
end
end

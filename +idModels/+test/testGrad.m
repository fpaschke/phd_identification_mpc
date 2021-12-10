function testGrad(fun,varargin)
%TESTGRAD Plot function output and the difference between nummerical and
%analytical Gradient.
% 
% testGrad(fun,x1,x2,...,xn,varargin)
% fun [function handle]:    function f that outputs its value and its gradient
%                           wrt specified parametervector (at least two outputs) 
% x1,...,xn [unspecified]:  inputs to fun in corresponding order
% varargin: Name Value Pairs
%   'ix_p' [pos int]:   Inputindex of parametervector. If ix_p=1 then x1 will be
%                       treated as parameter. Def. 2
%   'ix_out' [pos int]: Outputindex of gradient. Def. 2
%   'delta' [pos double]: Stepsize

[varargin, ix_p, ix_out ,delta] = locExtPar(varargin,'ix_p',2,'ix_out',2,'delta',1e-8);

p = varargin{ix_p};
out = cell(1,ix_out);
[out{:}] = fun(varargin{:});
y = out{1};
dy = out{ix_out};
subplot(size(p,2)+1,1,1); plot(y); 
for np = 1:size(p,2)
    dp = zeros(1,size(p,2));
    dp(np) = delta;
    %dy_= (fun(u,p+dp,varargin{:}) - fun(u,p-dp,varargin{:}))/(2*delta);
    dy_= (fun(varargin{1:ix_p-1},p+dp,varargin{ix_p+1:end}) - fun(varargin{1:ix_p-1},p-dp,varargin{ix_p+1:end}))/(2*delta);
    subplot(size(p,2)+1,1,np+1); plot(dy(:,np)-dy_);
end
end

function [arg,varargout] = locExtPar(arg,varargin)
    varargout = cell(1,length(varargin)/2);
    for np = 1:2:length(varargin)
        ix = cellfun(@(x) isequal(x,varargin{np}),arg);
        if any(ix)
            varargout{(np+1)/2} = arg{find(ix)+1};
            arg(find(ix)+(0:1)) = [];
        else
            varargout{(np+1)/2} = varargin{np+1};
        end
    end
end
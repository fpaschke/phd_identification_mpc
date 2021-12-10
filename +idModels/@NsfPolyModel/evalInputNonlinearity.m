function [fuy, dfuy_dalpha] = evalInputNonlinearity(obj,u,yd)
%Calculates the Inputsignals to the linear Block of the system by evaluating the inputnonlinearities.
%
% [fuy, dfuy_dalpha] = evalInputNonlinearity(obj,u,yd)
% u [cell of Ns x nu double arrays]: The inputSignalvector into the inputnonlinearity
% yd [cell of Ns x ny double arrays]: The outputsignal that will be fed to f(u,y) if System has an Outputfeedback.
% fuy [cell of Ns x nu_lti double arrays]: The inputsignals to the linear block of the system
% dfuy_dalpha [cell of Ns x np double arrays]: The Gradients (np is the number of parameters of the corresponding inputnonlinearity).

nu = size(obj.Nk,2); cellin = 1;
if ~iscell(u); u = {u}; cellin = 0; end
if nargin==3 && ~iscell(yd); yd = {yd}; end

Nsets = length(u);
if nargout>1; dfuy_dalpha = cell(Nsets,nu); end
fuy = cell(length(u),1);

f = obj.InputNonlinearity;
f_nouts = obj.F_nargout;
for ni = 1:nu
    pi = obj.getFiParams(ni);
    for ns = 1:Nsets
        out = cell(f_nouts(ni),1);
        fi = f(ni);
        if obj.HasInputNonlinearity(ni)  %Hammer-Sys
            if obj.HasOutputFeedback(ni) %NL-Feedback Sys
                [out{:}] = fi.fun(u{ns}(:,fi.input_idx),pi,yd{ns}(:,fi.output_idx));
            else %Hammerstein sys
                [out{:}] = fi.fun(u{ns}(:,fi.input_idx),pi);
            end
        elseif ~isempty(fi.input_idx) 
            [out{:}] = u{ns}(:,fi.input_idx);
        else
            error('Something is wrong with obj!');
        end
        if f_nouts(ni) == 1
            fuy{ns}(:,ni) = out{1};
        else
            [fuy{ns}(:,ni),dfuy_dalpha{ns,ni}] = out{1:2};
        end
    end
end

if ~cellin
    fuy = fuy{1};
end
end
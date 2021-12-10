function [x_hat,y_hat,e] = observe(obj,y,u,x0)
%OBSERVE Estimate states of NsfPolyModel using IO-Data.
% 
%[x_hat,y_hat,e] = observe(obj,y,u,x0)
% y [(cell of) N x ny double]: The Measured Outputs of the Model. Each column represents one output. y = [y(1); ...; y(N)]
% u [(cell of) N x nu double]: The Measured Inputs. Each column represents one input. u = [u(1); ...; u(N)]
% x0 [nx x Nset double]: Initial conditions at sampling instant ("index") 1 where Nset is the number of datasets.
% y_hat [(cell of) N x ny double]: The Observed Output of the Model. y_hat = [y_hat(1); ...; y_hat(N)]
% x_hat [(cell of) N x nx double]: The Observed state of the Model. y_hat = [x_hat(1)=X0; ...; x_hat(N)]
% e [(cell of) N x nx double]: The Obervation Error e: e = y - y_hat.

assert(~isempty(obj.Ss),'System isnt fully parameterized which means that any of the polynomials A,B,... or E contains NaNs!');
fu = obj.evalInputNonlinearity(u,y);
[y_hat, x_hat, e] = idModels.alg.lti.observe(obj.Ss,y,fu,x0,obj.OutputOffset); %xobs(1|0) = X0   ... %yobs(1)=ypred(1) ...
end


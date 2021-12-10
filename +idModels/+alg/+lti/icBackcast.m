function [e_hat,z0,w_hat] = icBackcast(C,w)
%Computes initial conditions for ARMA(X) process Ce = Ay - Bu = w of a general Hammerstein ARMA(X) Model where the input 
%w is assumed to be the sequence w = Ay - Bf(u,alpha). 
%
% [e_hat,z0,w_hat] = icBackcast(C,w)
% C [double vector]: Coffizients of monic polynomial C = 1 + c1 q^-1 + ... + cnc q^-nc 
% w [cell of double vectors or double vector]:  The sequence(s) w. Compute be MA-Filtering w = Ay-Bf(u,alpha). If w is a cell
%                                               the Alg. assumes Nset independet Datasets for which the ICs will be comp.
% e_hat [nc x Nset double matrix]: Init. Cond. e(ts-1),...,e(ts-nc) for each of the NSet sequences w, where ts max(na,nb)
% z0 [nc x Nset double matrix]: Init. Cond. for matlab filter implementation for each of the NSet sequences w. 
% w_hat [nc x Nset double matrix]: Backcasted values of w

if ~iscell(w) w = {w}; end
assert(C(1)==1,'C needs to be monic!');

Nset = length(w);
nc = length(C)-1;

% Comp necessary matrices
Hc = hankel(C(2:end));
Hh = hankel(C(end-1:-1:1));
Hh = Hh(:,end:-1:1);

for k = 1:Nset
    % 1. Step: Backwardfiltering of sequence w = Ce -> e = 1/C*w
    [eb_hat, ~] = filter(1,C,w{k}(end:-1:1)); %[eb_hat(N) eb_hat(N-1)... eb_hat(ts)]
%     %TEST
%     eb_hat2(length(w{k})+[1:nc],1) = 0;
%     for l = length(w{k}):-1:1
%         eb_hat2(l) = -C(2:end)*eb_hat2(l+1:l+nc) + w{k}(l);
%     end
%     %TEST END
    
    % 2.Step: Backforecast w(k) for k<ts 
    w_hat(:,k) = Hc*eb_hat(end:-1:end-nc+1); %w_hat(ts-1) ... w_hat(ts-nc)

    % 3. Step: Backforecast e(k) for k<ts 
    e_hat(:,k) = Hh\w_hat(:,k); % e_hat(ts-1) ... e_hat(ts-nc)
    if nargout > 1
        z0(:,k) = filtic(1,C,e_hat(:,k));
    end
end
end

%% Test
% C = conv([1 -.85+.3i],[1 -.85-.3i]);
% nc = length(C)-1;
% Nexp = 100;
% Nsamp = 50;
% 
% e0 = rand(nc,Nexp);
% e = rand(Nsamp,Nexp);
% w = filter(C,1,e,e0);
% [~,ic] = idModels.alg.lti.icBackcast(C,mat2cell(w,Nsamp,ones(1,Nexp)));
% e_rec = filter(1,C,w) - e;
% e_opt = filter(1,C,w,ic) - e;
% plot([var(e_rec,[],2) var(e_opt,[],2)]); 
% legend({'AW 0' 'AW - Backcast'}); ylabel('\sigma^2(k)'); title('\sigma^2(k) =var(e(k) - e_{rec}(k))'); xlabel('k')
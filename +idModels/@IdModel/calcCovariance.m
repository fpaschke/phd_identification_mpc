function Lambda = calcCovariance(e,Npar,SubMean)
%CALCCOVARIANCE Calculate Covariance Matrix of e. NaN's will be ignored.
%
% Lambda = calcCovariance(e,Npar)
% e [double Matrix or Cell of double Matrices]: Signals with samples row wise.
% Npar [opt. pos int]: If supplied, then the corrected estimate of Lamda will be used. See Söderström, Stoica Chap. 3. Def. 0
% SubMean [opt. log.]. If 1 the mean of e will be subtracted.

if nargin == 1 || isempty(Npar); Npar = 0; end
assert(mod(Npar,1)==0,'Npar needs to be a pos. integer!')

if nargin < 3 || isempty(SubMean); SubMean = false; end
assert(SubMean==0 || SubMean==1,'Npar needs to be a pos. integer!')

if iscell(e)
    e = cell2mat(e);
end

if SubMean
    e = e - mean(e,1,'omitnan');
end

ix_dat = all(~isnan(e),2);
Ns_red = sum(ix_dat,1);

% Cv2 = e(ix_dat,:)'*e(ix_dat,:); % irgendwie nicht das gleiche wie Cv???
e = e(ix_dat,:);
Cv = e'*e;
if Ns_red-Npar <=0 
    warning('Number of samples is to low! -> Using unbiased estimate of covariance!');
    Npar = 0;
end
Lambda = Cv/(Ns_red-Npar);
end
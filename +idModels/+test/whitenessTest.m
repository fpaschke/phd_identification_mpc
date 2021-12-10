function [pval,H0_acc,t] = whitenessTest(x,m,alpha,combined)
%Check sequence(s) x for whiteness. The test uses Autocorrelation series Rxx(k) of x to compute teststatistic 
%T = N/Rxx(0)^2*(R(1)^2 ... R(m)^2). The Teststatistic should then be asymptotically Chi-squared distributed with m degrees  
%of freedom in case of white noise.
%Nullhyp. H0: "x is white noise".
%
% [pval,H0_acc,t] = whitenessTest(x,m,alpha,combined)
% x [double vector or cell of double vectors]: The sequence(s) wich should be tested.
% m [pos. int. scalar]: Number of lags to compute teststatistic.
% alpha [opt. double scalar 0...1]:  Significance niveau (lower bound of pvalue). Determines the area of the distribution of the Teststatistic which is
%                               regarded as the critical ("over-random") area. 0.05 eg. means that the uppermost 5% of the 
%                               distribution of the Teststatistic are regarded as critical. If the the p-value of the
%                               computed teststatistic is in the critical area the test fails and returns 0, which means H0
%                               will be rejected. Typical values are 0.05, 0.01 and 0.001. Default: alpha == .01
% combined [opt 0/1]:   If multiple datasets are used and combined == true then the algorithm puts all datasets together
%                       and determines whiteness of combined sequence. Thus only one p value and one testresult will be
%                       computed.
% H0_acc [logical]: true indicates that there is no reason to reject H0 wrt the given significance alpha
% pval [double 0 ... 1]: the pvalue: the likelihood of observing such or a more extreme result under H0 i.e P(T>=t|H0). If 
%                        p is 0 H0 should be rejected. If p is 1 then there is no reason to reject H0.
% t [double scalar]: value of teststatistic
%                               
% EXAMPLES:
% X = randn(1,1000);
% Y = filter(1,[1 -0.5],X);
% 
% [p1,f1] = bModels.test.whitenesstest(X,50,0.01)
% [p2,f2] = bModels.test.whitenesstest(Y,50,0.01)
% [p3,f3] = bModels.test.whitenesstest({Y X},50,0.01,1)
% [p4,f4] = bModels.test.whitenesstest({Y X},50,0.01)

if nargin<3 || isempty(alpha)
    alpha = .01;
end
if nargin<4 || isempty(combined)
    combined = false;
end
assert(combined==0 || combined==1,'The parameter "combined" need to be either false or true!');
if combined == true
   x = cell2mat(x); 
end

if ~iscell(x)
    x = {x};
end
Nset = length(x);

H0_acc = false(Nset,1);
pval = NaN(Nset,1); t = NaN(Nset,1);
for ns = 1:Nset
    assert(any(size(x{ns})==1),'The sequences need to be a vectors');
    N = length(x{ns});
    if m > length(x{ns}) - 1
        mm = length(x{ns}) - 1;
    else
        mm = m;
    end
	[XC,lags] = xcorr(x{ns}(:),mm,'biased');  
    l0 = find(lags == 0);
    %plot(lags,XC/XC(l0)); hold on; plot(mm*[-1 1],3/sqrt(N)*[1 1; -1 -1]','--r')
    [pval(ns), H0_acc(ns),t(ns)] = calcAndDecide(XC(l0+1:end),XC(l0));
end


%% INLINE
function [p,r,t] = calcAndDecide(xc,sigma2)
    t = N/sigma2^2*sum(xc.^2);
    try %statistic toolbox available
        p = 1 - chi2cdf(t,length(xc));
    catch 
        p = 1 - gammainc(t/2,length(xc)/2);
    end
        
    if p<alpha
        r=false;
    else 
        r=true;
    end
end
end
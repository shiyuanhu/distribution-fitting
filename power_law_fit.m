function [mu,C,AIC] = power_law_fit(x,a,b,mu0)
% Fit x to a power-law distribution using a maximum likelihood 
% estimate. Assume x follows phi = @(x) c*x^{-mu}, where c is the 
% normalization constant. x is distributed between a and b. 
% a is the low limit, typically taken to be the min(x)
%
% Input: 
%       x: 1D array of the random variable to be fitted
%       a, b: Minimum and maximum values of the random variable.
%             If b = inf, then there is no upper bound. 
%       mu0: initial guess of mu
% 
% Return: 
%       mu: the best estimate parameter
%       C: a function handle to compute the complement of the 
%          cumulative distribution
%       AIC: a second-order Akaike's Information Criterion that 
%            maybe later used for model selection
%
% Written by Shiyuan Hu <shiyuan.hu@nyu.edu>, Dec. 17, 2019
%
N = length(x);
%
% form negative log likelihood
minus_L = @(mu)(-lnlike_power(mu,x,a,b));
%
% minimize the negative log likelihood
[mu,max_L] = fminsearch(minus_L,mu0);
%
% form C
C = @(y) 1-(a^(1-mu)-y.^(1-mu))./(a^(1-mu)-b^(1-mu));
%
% compute AIC
K = 2;
AICs_power = 2*max_L+2*K+2*K*(K+1)/(N-K-1);

    function L = lnlike_power(mu,x,a,b)
        % log likelihood function
        N = length(x);
        L = N*log(mu-1)-N*log(a^(1-mu)-b^(1-mu))-mu*sum(log(x));
    end

end

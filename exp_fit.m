function [theta,C,AIC] = exp_fit(x,a,b,theta0)
% Fit x to an exponential distribution using a maximum likelihood
% estimate. Assume x is identically and independently from 
% phi = @(x) c*exp(-theta x), where c is the normalization constant.
% Assume x is distributed between a and b, then 
% c = theta/(exp(-a*theta)-exp(-b*theta)).
% 
% A distribution can be visualized without histogram by the 
% rank-frequency distribution: C(y) = 1-\int_a^y phi(x)dx, which is the 
% complement of the culmulative distribution.
%
% Input:
%       x: 1D array of the random variables to be fitted
%       a, b: Minimum and maximum values of the random variable.
%             Generally, a can be taken as min(x) and b can be take as
%             max(x), or inf when there is no upper bound.
%       theta0: initial guess of theta 
% Return:
%       theta: the best estimate parameter
%       C: a function handle to compute C(x)
%       AIC: a second-order Akaike’s Information Criterion that maybe later
%            used for model selction
%  
% Written by Shiyuan Hu <shiyuan.hu@nyu.edu>, Dec. 17, 2019
%
N = length(x);
%
% form negative log likelihood
minus_L_exp = @(theta)(-lnlike_exp(theta,x)); 
%
% minimize the negative log likelihood
[theta, max_L] = fminsearch(minus_L_exp,theta0);
% 
% form C
c = theta/(exp(-a*theta)-exp(-b*theta));
C = @(y)1-c/theta*(exp(-a*theta)-exp(-y*theta));
%
% compute AIC
K = 2;
AIC = 2*max_L+2*K+2*K*(K+1)/(N-K-1);

    function L = lnlike_exp(theta,x)
        c = theta/(exp(-a*theta)-exp(-b*theta));
        L = N*log(c)-theta*sum(x);
    end
end

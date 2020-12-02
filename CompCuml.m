function [C, x] = CompCuml(x)
% Compute the Complement of culmulative distribution of 1D array of random 
% variable.
% 
% Input: 
%       x: 1D array of the random variables
% Return:
%       C: the Complement of culmulative distribution
%       x: sorted x
%
    x = sort(x,'ascend');
    C = (length(x):-1:1)./length(x);
end
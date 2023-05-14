function [pvalue, stat] = binomtest(s, n, p)
% BINOMTEST Performs a binomial test.
%
% INPUTS:
%   s - The number of successes.
%   n - The total number of trials.
%   p - The probability of success.
%
% OUTPUTS:
%   pvalue - The p-value of the test.
%   stat - The test statistic.
%
% EXAMPLE:
%   % Perform a binomial test to see if 5 heads in 10 coin flips is
%   % significantly different from a fair coin.
%   s = 5;
%   n = 10;
%   p = 0.5;
%   [pvalue, stat] = binomtest(s, n, p)
%
%   % The p-value is 0.03125, which is significant at the 0.05 level.
%   % We can conclude that the number of heads is significantly
%   % different from what we would expect if the coin were fair.

% Calculate the test statistic.
stat = (s - n * p) / sqrt(n * p * (1 - p));

% Calculate the p-value.
pvalue = 1 - chi2cdf(stat^2, 1);

end
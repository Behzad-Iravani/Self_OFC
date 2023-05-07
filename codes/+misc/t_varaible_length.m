% -*- UFT-8 -*-
function [t, all, n, s, ci] = t_varaible_length(T, lock_type,s)
%   t_varaible_length calculates the statistics for table T
%   Input:
%          -T:          a table containing data
%          -lock_type:  a string or character vector determines the type of lock (stim or resp)
%          -s:          a scalar defines the amount of smooting
%   Output:
%          -t:          a vector of t-scores for each time point in the data.
%          -all:        a matrix of smoothed and zero-padded data for each trial in the table.
%          -n:          a vector of trial counts for each time point in the data.
%          -s:          a vector of variance estimates for each time point in the data.
%          -ci:         a 2xM matrix containing the lower and upper confidence interval bounds for each time point, where M is the length of the maximum trial.
%   
%   Copyright (C)  Behzad Iravani, department of neurology and neurological
%   sciences, Stanford University. May 2023
%   Author: Behzad Iravani
%   behzadiravani@gmail.com
%   Contact: behzadiravani@gmail.com
%   Date: 05/02/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments % checks the arguments
    T table % a table of data
    lock_type string  % the type of lock ('stim' or 'resp')
    s double {mustBeGreaterThan(s, 0)} % smoothing parameter (must be greater than zero)
end

if nargin < 3 % checks if s is provided or not
    s = .1; % default smoothing parameter value
    fprintf('The smoothing parameter is set to %1.2f\n', s); % display a message indicating the default value used
end

T = stat.average_over_task(T); % average over the task

% determining the maximum length of trials
mxt = max(cellfun(@length,T.avg)); % mxt: the maximum number of time points in all trials
all = zeros(height(T), mxt); % all: a matrix to store the averaged data for all trials, with the number of rows equal to the number of trials and the number of columns equal to mxt
tmp = zeros(1, mxt); % tmp: a temporary vector to store the sum of all data points for each time point
n   = zeros(1, mxt); % n: a vector to store the number of data points for each time point
%------------------------------------
baseline = []; % initializing the baseline
if strcmp(lock_type, 'resp') % checking the lock_type
    % --- response lock 
    for h=1:height(T) % loops over the rows of T
        tmpa = zeros(1, mxt); % initialize a temporary vector with zeros of length mxt
        % Shift the current row's trial data to the end of the tmpa vector,
        % leaving zeros in the beginning. This is done to align all trials
        % to the end of the vector for averaging
        tmpa(end-length(T.avg{h})+1:end) = T.avg{h};

        % Calculate the baseline value for the current trial, which is the
        % mean of the trial's data after time = 0, for response lock and
        % before time = 0, for stimulus lock
        baseline(h,:) = mean(T.avg{h}(T.time{h}>0)); % the baseline for the response lock is after zero
        
        % Smooth the trial data using a moving average filter of window size
        % specified by the input parameter s. considering the fsample 512 for
        % defalut s value the amound for smoothig is .1 x 1000 /512 ~ 200ms smoothing 
        all(h,:) = movmean(tmpa ,floor(s*1e3));

        % Add the current trial's data to the temporary vector tmp.
        % tmp is used to calculate the mean of all trials.
        tmp = tmp + tmpa ;

        % Create a temporary vector of ones with length mxt and shift it to
        % the end of the vector to align with the current trial. This vector
        % is used to keep track of the number of valid trials for each time point.
        tmpn = zeros(1, mxt);
        tmpn(end-length(T.avg{h})+1:end) = 1;
        n = n + tmpn;% add the tmpn to the n vector
    end % for loop over T's row 

elseif strcmp(lock_type, 'stim')
    % --- stimulus lock 
    for h=1:height(T)% loops over the rows of T
        tmpa = zeros(1, mxt); % initialize a temporary vector with zeros of length mxt
        tmpa(1:length(T.avg{h})) = T.avg{h};
        baseline(h,:) = mean(T.avg{h}(T.time{h}<0));
        all(h,:) = movmean(tmpa ,floor(s*1e3));
        tmp = tmp + tmpa;


        tmpn = zeros(1, mxt);
        tmpn(1:length(T.avg{h})) = 1;
        n = n + tmpn;
    end % for h
end % if lock type
% compute statistics 
if n>1 % checks if more than one datapoint exists to perform statisitcs
m = tmp./n;
s  = sum((all-repmat(m, size(all,1),1)).^2)./(n-1);

ci(1,:) = m - 1.96*sqrt(s);
ci(2,:) = m + 1.96*sqrt(s);
t = sqrt(n).*(m-mean(baseline))./sqrt(s);
else
m = tmp./n;
s  = sum((all-repmat(m, size(all,1),1)).^2)./(n-1);

ci(1,:) = m ;
ci(2,:) = m ;

t = m;

end
end % t_varaible_length
% $END






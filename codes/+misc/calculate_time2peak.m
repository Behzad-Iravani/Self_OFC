% -*- UTF-8 -*-
function time2peak = calculate_time2peak(EP,SJ)
%   calculate_time2peak estimates the peak latecy for HFB enevelope across the 
%   anatomical sites. 
%   Input:
%      - EP: a table containing data for self-episodic condition 
%      - SJ: a table containing data for self-judgment condotion 
%   Output:
%      - time2peak: a matrix containigs the time (% RT) for the prominent
%      peak of HFB for every subjects
%
%   calculate_time2peak is part of the scripts that recreates the plots 
%   that were reported in "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
%   

% -----------------------------------------------------------


Tind      = [EP;SJ]; % concatenates the EP and SJ table
if iscategorical(Tind.task)% check if the task cloumn is categorical 
    Tind.task = cellstr(Tind.task); % converts back to string 
end % if categorical 
s = .1; % smoothing factor

isubj = 0; % set the subject counter for the for loop  
for subj = unique(EP.subj)' % loop over all the unique subjects in the data 
    isubj = isubj +1; % increament the isubj counter 
    for JP = ["OFC","MPFC"] % loop over the anatomical sites 
        % find the index for a given subject, given anatomical site and
        % self-referential tasks (self-episodic + self-judgment)
        index = cellfun(@(x) any(strcmp(x, JP)), Tind.JPAnatomy)...
            & strcmp(Tind.subj, subj{:})...
            & (strcmp(Tind.task, "EP") | strcmp(Tind.task, "SJ"));
        % computing the t-statistic 
        [time_series, all.(JP), nn.(JP), ss.(JP), ci_.(JP)] = misc.t_varaible_length(Tind(index,:), 'stim',s);
        % z-score the time-course 
        time_series  = zscore(time_series);
        % find the peak for the smoothed data using matlab findpeaks algorithm  
        [~,locs,~,p] = findpeaks(movmean(time_series,floor(s*1e3)));
        % find the peak with the maximum prominence 
        [~, locs_p]  = max(p);
        % create a time axis as percentage of reaction time 
        time         = linspace(0,100, length(time_series));
        % storing time ro peak for the give site
        time2peak.(JP)(isubj) = time(locs(locs_p));
    end % for sites JP
end % subjs
% $END

%-*- UTF8 -*-
function T_ = average_over_sessions(T)
%   average_over_session calculates the mean of Tval and Pval over the
%   sesions for each unique combination of channel, subject and tasks. 
%   Input:
%      - T: table contains data including subject name, channel label,
%      channel' MNI coordinate and etc. see stat_report for more
%      information. 
%   average_over_sessions is part of the scripts that recreates the plots and
%   that was reported in "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
%   
%   Copyright (C)  Behzad Iravani, department of neurology and neurological
%   sciences, Stanford University. May 2023 
%
%   Author: Behzad Iravani
%   behzadiravani@gmail.com
%   Contact: behzadiravani@gmail.com
%   Date: 05/02/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscategorical(T.subj) % chack if subj is categorical 
   T.subj = cellstr(T.subj); % converts it to cell string 
end
if iscategorical(T.task) % chack if task is categorical 
   T.task = cellstr(T.task); % converts it to cell string 
end

[~, ~, IC] = unique(strcat(T.chan, ':', T.subj, ':', T.task)); % find a unique combinations of channel, subjects and tasks
T_ = table();  % sets up an empty output table T_
warning off
% find the index of columns to be included in the output table 
cindex = 0;  % set the cindex that is used as a counter in the following for loop  
for columnname = ["subj", "chan", "task", "X", "Y", "Z", "dof", "responseTr", "Loc_Tval", "Loc_Pval", "avg", "time", "RT", "varRT", "JPAnatomy", "Density", "BDI", "BDA", "absX"]% loop over the column names
    if ~isempty(find(strcmp(T.Properties.VariableNames, columnname))) % if the column exists in the input table store the index of column in the variable index 
        cindex =  cindex +1;  % increament this for loop counter 
        index(cindex) = find(strcmp(T.Properties.VariableNames, columnname)); % store the index of the current column
    end
end
% -------
for ri = 1:max(IC) % loop over all the unique rows in the input table 
    T_(ri,1:length(index)) = T(find(IC == ri,1,'first'),[index]); % find the first row of the given non-unique rows
    T_.Tval(ri) = nanmean(T.Tval(IC == ri)); % average t-value for baseline comparison 
    T_.Pval(ri) = nanmean(T.Pval(IC == ri)); % average p-value for baseline comparison 

end
warning on% turn back the warning on
% adding the name of the new columns
T_.Properties.VariableNames(1:(length(index)+2)) = [T.Properties.VariableNames([index]), {'Tval'}, {'Pval'}];
end % average_over_sessions
% $ END
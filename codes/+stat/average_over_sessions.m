%-*- UTF8 -*-
function T_ = average_over_sessions(T)
% average_over_session calculates the mean of Tval and Pval over the
% sesions for each unique combination of channel, subject and tasks. 
% Input:
%      - T: table contains data including subject naem, channel label,
%      channel' MNI coordinate and etc. see stat_report for more
%      information. 
%   average_over_sessions is part of the scripts that recreates the plots and
%   that was reported in "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
%   Author: Behzad Iravani
%   behzadiravani@gmail.com
%   Contact: behzadiravani@gmail.com
%   Date: 05/02/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, ~, IC] = unique(strcat(T.chan, ':', T.subj, ':', T.task)); % find a unique combinations of channel, subjects and tasks
T_ = table();  % sets up an empty output table T_
warning off
% find the index of columns to be included in the output table 
cindex = 0;
for columnname = ["subj", "chan", "task", "X", "Y", "Z", "dof", "responseTr", "Loc_Tval", "Loc_Pval", "avg", "time", "RT", "varRT", "JPAnatomy", "Density", "BDI", "BDA"]
    if ~isempty(find(strcmp(T.Properties.VariableNames, columnname)))
        cindex =  cindex +1 ;
        index(cindex) = find(strcmp(T.Properties.VariableNames, columnname));
    end
end
% -------
for ri = 1:max(IC)
    T_(ri,1:length(index)) = T(find(IC == ri,1,'first'),[index]);
    T_.Tval(ri) = nanmean(T.Tval(IC == ri)); % average t-value for baseline comparison 
    T_.Pval(ri) = nanmean(T.Pval(IC == ri)); % average p-value for baseline comparison 

end
warning on
% adds the columns' name to the output table  
T_.Properties.VariableNames(1:(length(index)+2)) = [T.Properties.VariableNames([index]), {'Tval'}, {'Pval'}];
end
% $ END
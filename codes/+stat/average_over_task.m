%-*- UTF8 -*-
function T_ = average_over_task(T)
% average_over_task calculates the mean of Tval and Pval over the
% task for each unique combination of channel, and subjects.
% Input:
%      - T: table contains data including subject name, channel label,
%      channel' MNI coordinate and etc. see stat_report for more
%      information.
%   average_over_sessions is part of the scripts that recreates the plots and
%   that was reported in "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
%
%   Copyright (C)  Behzad Iravani, department of neurology and neurological
%   sciences, Stanford University. May 2023.
%
%   Author: Behzad Iravani
%   behzadiravani@gmail.com
%   Contact: behzadiravani@gmail.com
%   Date: 05/02/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if ses exists in the table
assert(~ismember({'ses'}, T.Properties.VariableNames), 'Data should be first averaged over session by calling average_over_session')

[~, ~, IC] = unique(strcat(T.chan, ':', T.subj)); % find a unique combinations of channel and subjects

T_ = table(); % sets up an empty output table T_
warning off
% find the index of columns to be included in the output table
cindex = 0;  % set the cindex that is used as a counter in the following for loop
for columnname = ["subj", "chan", "X", "Y", "Z", "dof", "time", "JPAnatomy", "Density", "BDI", "BDA", "absX"]
    if ~isempty(find(strcmp(T.Properties.VariableNames, columnname)))
        cindex =  cindex +1 ;
        index(cindex) = find(strcmp(T.Properties.VariableNames, columnname));
    end
end
% -------
df_c = 4; % column to be added in the end
additional_column_name ={};
for ri = 1:max(IC)
    T_(ri,1:length(index)) = T(find(IC == ri,1,'first'),[index]);
    T_.Tval(ri)            = mean(T.Tval(IC == ri)); % average t-value for baseline comparison
    T_.Pval(ri)            = mean(T.Pval(IC == ri)); % average p-value for baseline comparison


    T_.RT(ri)              = mean(T.RT(IC == ri)); % average reaction time (RT)
    T_.varRT(ri)           = mean(T.varRT(IC == ri)); % average variance of reaction time (RT)
    if ismember({'responseTr'}, T.Properties.VariableNames)
        T_.responseTr(ri)      = mean(T.responseTr(IC == ri));
        if ri ==1
            df_c = df_c +1;
            additional_column_name{end+1} =  {'responseTr'};
        end
    end
    if ismember({'Loc_Tval'}, T.Properties.VariableNames)
        T_.Loc_Tval (ri)       = mean(T.Loc_Tval(IC == ri)); % average p-value for localization (Self>Math)
        if ri ==1
            df_c = df_c +1;
            additional_column_name{end+1} =  {'Loc_Tval'};
        end
    end
    if ismember({'Loc_Pval'}, T.Properties.VariableNames)
        T_.Loc_Pval(ri)        = mean(T.Loc_Pval(IC == ri));  % average t-value for localization (Self>Math)
        if ri ==1
            df_c = df_c +1;
            additional_column_name{end+1} =  {'Loc_Pval'};
        end
    end
    if  ismember('avg', T.Properties.VariableNames) % if avg exists average it over the tasks
        if sum(IC == ri)>1
            T_.avg{ri}           = mean(cat(1,T.avg{IC == ri}));
        else
            T_.avg{ri}           = cat(1,T.avg{IC == ri});
        end
        if ri ==1
            df_c = df_c +1;
            additional_column_name{end+1} =  {'avg'};
        end
    end

end
warning on
% adds the columns' name to the output table
T_.Properties.VariableNames(1:(length(index) + df_c)) = [T.Properties.VariableNames([index]),...
    {'Tval'}, {'Pval'}, {'RT'}, {'varRT'}, additional_column_name{:}];


end
% $ END
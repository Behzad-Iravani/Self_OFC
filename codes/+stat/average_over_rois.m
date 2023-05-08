% -*- UTF-8 -*-
function T_ = average_over_rois(T)
%   average_over_rois calculates the mean of Tval and Pval over the
%   sesions for each unique ROIs.
%   Input:
%      - T: table contains data including subject name, channel label,
%      channels' MNI coordinate and etc. see stat_report for more
%      information.
%   average_over_rois is part of the scripts that recreates the plots
%   that were reported in "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
%
%   Copyright (C)  Behzad Iravani, department of neurology and neurological
%   sciences, Stanford University. May 2023.
%
%   Author: Behzad Iravani
%   behzadiravani@gmail.com
%   Contact: behzadiravani@gmail.com
%   Date: 05/06/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[C, ~, IC] = unique(T.JPAnatomy); % finds the indicies for unique ROIs in the data
T_         = table(); % initiate the output tabel T_
warning off % turns the warning off

cindex = 0; % set the cindex that is used as a counter in the following for loop
for columnname = ["subj", "chan", "task", "dof", "responseTr", "time", "RT", "varRT", "Density", "BDI", "BDA"] % loop over the column names
    if ~isempty(find(strcmp(T.Properties.VariableNames, columnname))) % if the column exists in the input table store the index of column in the variable index
        cindex =  cindex +1; % increament this for loop counter
        index(cindex) = find(strcmp(T.Properties.VariableNames, columnname)); % store the index of the current column
    end % if
end % for
df_c = 9; % column to be added in the end
additional_column_name ={};
for ri = 1:max(IC) % loop over all the unique rows in the input table

    T_(ri,1:length(index)) = T(find(IC == ri ,1,"first"),[index]); % find the first row of the given non-unique rows
    T_.JPAnatomy(ri)       = C(ri); % storing the name of region
    % calculates the descriptive statistics:
    T_.Tval(ri)            = mean(T.Tval(IC == ri));
    T_.TvalMed(ri)         = median(T.Tval(IC == ri)); % added this on 5/3/2023
    T_.Pval(ri)            = mean(T.Pval(IC == ri));
    % calucltes the localizing statistics:
    T_.Loc_Tval(ri)            = mean(T.Loc_Tval(IC == ri));
    T_.Loc_Pval(ri)            = mean(T.Loc_Pval(IC == ri));
   % calculates the median of all the electrode's coordinates that contributed to statistics
    T_.X(ri) = median(T.X(IC == ri)); % added this on 5/3/2023
    T_.Y(ri) = median(T.Y(IC == ri)); % added this on 5/3/2023
    T_.Z(ri) = median(T.Z(IC == ri)); % added this on 5/3/2023
    if ismember({'Tval_pred'}, T.Properties.VariableNames)
        T_.Tval_pred(ri)            = mean(T.Tval_pred(IC == ri));
        if ri ==1
        df_c = df_c +1;
        additional_column_name{end+1} =  {'Tval_pred'};
        end
    end
    if ismember({'Tval_predL'}, T.Properties.VariableNames)
        T_.Tval_predL(ri)            = mean(T.Tval_predL(IC == ri));
        if ri ==1
        df_c = df_c +1;
        additional_column_name{end+1} ={'Tval_predL'};
        end
    end
    if ismember({'Tval_predH'}, T.Properties.VariableNames)
        T_.Tval_predH(ri)            = mean(T.Tval_predH(IC == ri));
        if ri ==1
        df_c = df_c +1;
        additional_column_name{end+1} ={'Tval_predH'};
        end
    end
  
    % check if avg exists in the data
    if ismember({'avg'}, T.Properties.VariableNames)
        if sum(IC == ri)>1 % checks if more than one data points exist
            T_.avg{ri}           = mean(cat(1,T.avg{IC == ri}));
        else
            T_.avg{ri}           = cat(1,T.avg{IC == ri});
        end
        if ri ==1
        df_c = df_c +1;
        additional_column_name{end+1} = {'avg'};
        end

    end
 
end % for ri
warning on % turn back the warning on
% adding the name of the new columns
T_.Properties.VariableNames(1:(length(index) + df_c)) = [T.Properties.VariableNames([index]),...
    {'JPAnatomy'}, {'Tval'}, {'TvalMed'}, {'Pval'}, {'Loc_Tval'}, {'Loc_Pval'}, {'X'}, {'Y'}, {'Z'},  additional_column_name{:}];


end % average_over_rois
% $END
function T_ = average_over_sessions(T)


[~, ~, IC] = unique(strcat(T.subj));
T_ = table();
warning off
cindex = 0;
for columnname = ["subj", "dof", "responseTr", "RT", "varRT", "JPAnatomy", "Density", "BDI", "BDA"]
    if ~isempty(find(strcmp(T.Properties.VariableNames, columnname)))
        cindex =  cindex +1 ;
        index(cindex) = find(strcmp(T.Properties.VariableNames, columnname));
    end
end

for ri = 1:max(IC)
    T_(ri,1:length(index)) = T(find(IC == ri,1,'first'),[index]);
    T_.Tval(ri)    = nanmean(T.Tval(IC == ri));
    T_.TvalMed(ri) = nanmedian(T.Tval(IC == ri));
    T_.Pval(ri)    = nanmean(T.Pval(IC == ri));



    T_.X(ri) = median(T.X(IC == ri)); % added this ub 5/3/2023
    T_.Y(ri) = median(T.Y(IC == ri)); % added this ub 5/3/2023
    T_.Z(ri) = median(T.Z(IC == ri)); % added this ub 5/3/2023

    T_.CvmPFC(ri) = sum(strcmp(T.JPAnatomy(IC == ri), 'MPFC'));

end
warning on
T_.Properties.VariableNames(1:(length(index)+7)) = [T.Properties.VariableNames([index]),...
    {'Tval'}, {'TvalMed'} {'Pval'}, {'X'}, {'Y'}, {'Z'}, {'nMPFC'}];

end
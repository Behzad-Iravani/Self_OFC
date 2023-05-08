%-*-UTF-8 -*-
% Plot Fig 1E showing the difference between OFC and vmPFC activity for two groups of subjects
function plot_Fig1E(EP,SJ, col)
%   plot_Fig1E genertate the bar graphs that compare the EP(SE) and SJ HFB in the same brain
%   Input:
%         - EP:     a table containing the data for self-episodic 
%         - SJ:     a table containing the data for self-judgment
%         - col:    a 2X3 color matrix 
%
%   plot_Fig1E is part of the scripts that recreates the plots
%   that was reported in "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
%
%   Copyright (C)  Behzad Iravani, department of neurology and neurological
%   sciences, Stanford University. May 2023
%
%   Author: Behzad Iravani
%   behzadiravani@gmail.com
%   Contact: behzadiravani@gmail.com
%   Date: 05/08/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a new figure
figure
% Hold the plot
hold on
% Initialize variables for counting the number of conditions and subjects
bc = 0;
nc = 0;
% Create two subplots, one for each group of subjects
ax(1) = subplot(121);
xlim([-1 3]) % set the x-axis limits for axis 1
ylim([-2.5 .5])% set the y-axis limits for axis 1
axis square
title('\rmvmPFC-OFC') % adding the title for axis 1

ax(2) = subplot(122);
xlim([-1 4]) % set the x-axis limits for axis 2
ylim([-.5 3]) % set the y-axis limits for axis 2
axis square

% find the unique subjects in the both tabels 
us = unique([EP.subj; SJ.subj]);
% get unique colors per subjects 
colsubs = lines(length(us));
M = containers.Map(us, 1:length(us));

result=struct();

icolor = 0;
% Loop through the two groups of subjects
for tsk = ["EP", "SJ"]
    % Increment the condition counter
    nc = nc +1;
    % Evaluate expressions to extract values and subject IDs for MPFC and OFC
    eval(strcat("valm = ",tsk, ".Tval_pred(", tsk,".JPAnatomy == ""MPFC"")"));
    eval(strcat("valo = ",tsk, ".Tval_pred(", tsk,".JPAnatomy == ""OFC"")"));
    eval(strcat("valml = ",tsk, ".Tval_predL(", tsk,".JPAnatomy == ""MPFC"")"));
    eval(strcat("valol = ",tsk, ".Tval_predL(", tsk,".JPAnatomy == ""OFC"")"));
    eval(strcat("valmh = ",tsk, ".Tval_predH(", tsk,".JPAnatomy == ""MPFC"")"));
    eval(strcat("valoh = ",tsk, ".Tval_predH(", tsk,".JPAnatomy == ""OFC"")"));
    eval(strcat("subjm = ",tsk, ".subj(", tsk,".JPAnatomy == ""MPFC"")"));
    eval(strcat("subjo = ",tsk, ".subj(", tsk,".JPAnatomy == ""OFC"")"));
    % Compute mean difference between OFC and MPFC activity for each subject group
    result.mean(nc) = mean(valo - valm);
    result.low(nc)  = mean(valol - valml);
    result.high(nc)  = mean(valoh - valmh);
    result.std(nc) = abs((result.high(nc) - result.low(nc)) /(2 * 1.96))*sqrt(length(valm));
    % Loop through the values for each subject and plot them
    for ival = 1:length(valm)
        icolor= icolor +1;
        % Plot bar chart showing mean difference between OFC and MPFC activity for each subject group
        cx = bc;
        rn = randn(1,1);
        axes(ax(1))
        hold on
        bar(cx, [mean(valo - valm)], .8, 'FaceColor', col(nc,:), 'LineWidth', .75)
        errorbar(cx, mean(valo - valm),...
            abs( mean(valo - valm)- mean(valol - valml)),...
            abs(mean(valo - valm)- mean(valoh - valmh)), ...
            'LineStyle', 'none', 'CapSize', 0, 'LineWidth', 1.25, 'Color','k')
        hold off
        axes(ax(2))
        hold on
        scatter(cx+rn*.1, valo(ival), 'MarkerFaceColor', colsubs(M(subjm{ival}),:), 'MarkerEdgeColor', colsubs(M(subjm{ival}),:))

        cx = cx + 1;
        scatter(cx+rn*.1, valm(ival), 'MarkerFaceColor', colsubs(M(subjo{ival}),:), 'MarkerEdgeColor', colsubs(M(subjo{ival}),:))

        h = line([cx-1+rn*.1, cx+rn*.1],[valo(ival) valm(ival)],...
            'color', colsubs(M(subjm{ival}),:), 'LineStyle', '--', 'LineWidth', 1.25);
        hold off
    end
    bc = bc +2;
end
set(ax(1),'XTick', [0, 2],...
    'XTickLabel', {'SE' 'SJ', 'OFC'},...
    'YTick', [-2.5, 0, .5], 'LineWidth', 1.5, 'FontName', 'Arial', 'FontSize', 16)

set(ax(2),'XTick', unique(sort([.5:2:3,(.5:2:3)-.5 ,(.5:2:3)+.5])),...
    'XTickLabel', {'OFC', 'SE', 'vmPFC', 'OFC', 'SJ', 'vmPFC'},...
    'YTick', [-.5, 0, 3], 'LineWidth', 1.5, 'FontName', 'Arial', 'FontSize', 16)

ax(1).TickLength(1) = .025;
ax(2).TickLength(1) = .025;
box off % remove the box around the plot
axes(ax(1)) % select axis 1
title('\rmOFC-vmPFC') % add title to the plot
% print the plot to file
print -dsvg  results\Fig1E.svg
print -dpng -r300  results\Fig1E.png
% write the result to json file
result.p = 2*(1-tcdf(diff(result.mean)./...
    ((diff(result.high) - diff(result.low))/(2*1.96)),length(valm)-1));
result.stdc = (diff(result.high) - diff(result.low))/(2*1.96)*sqrt(length(valm));
json_txt = jsonencode(result, "PrettyPrint",true);
fid = fopen('vmPFCOFC_same_brain_stimlock.json', 'w');
fprintf(fid, json_txt);
fclose(fid);


end
% $END
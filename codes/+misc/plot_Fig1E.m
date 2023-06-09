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
xlim([-2 4]) % set the x-axis limits for axis 1
ylim([-1.5 4.5])% set the y-axis limits for axis 1



ax(2) = subplot(122);
xlim([-1 4]) % set the x-axis limits for axis 2
ylim([-.5 3]) % set the y-axis limits for axis 2


% find the unique subjects in the both tabels
us = unique([EP.subj; SJ.subj]);
% get unique colors per subjects
colsubs = lines(length(us));
M = containers.Map(us, 1:length(us));

result=struct();
dic = containers.Map([-1,0,1],1:3); % a dictionary for converting numbers to index
color_lines = arrayfun(@(j) misc.hex2rgb(j),["#C60071", "#5E5E5E", "#0073A7"],'UniformOutput',false); % define color for lines
% concat the colors
color_lines = cat(1,color_lines{:});
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
    result.mean(nc) = mean(valm - valo);
    result.low(nc)  = mean(valml - valol);
    result.high(nc)  = mean(valmh - valoh);
    result.std(nc) = abs((result.high(nc) - result.low(nc)) /(2 * 1.96))*sqrt(length(valm));
    cx = [0,2];
    for ival = 1:length(result.mean)
    axes(ax(1))
        hold on
        bar(cx(ival), result.mean(ival), 1.75, 'FaceColor', col(ival,:), 'LineWidth', .75)
        errorbar(cx(ival), result.mean(ival),...
            abs(result.low(ival)),...
            abs(result.high(ival)), ...
            'LineStyle', 'none', 'CapSize', 0, 'LineWidth', 1.25, 'Color','k')
        hold off
    end
    % Loop through the values for each subject and plot them
    for ival = 1:length(valm)
        icolor= icolor +1;
        % Plot bar chart showing mean difference between OFC and MPFC activity for each subject group
        cx = bc;
        rn = randn(1,1);
        axes(ax(2))
        hold on
        scatter(cx+rn*.1, valo(ival), 90, ...
            'MarkerFaceColor', colsubs(M(subjm{ival}),:), 'MarkerFaceAlpha', .45,...
            'MarkerEdgeColor', colsubs(M(subjm{ival}),:))

        cx = cx + 1;
        scatter(cx+rn*.1, valm(ival), 90,...
            'MarkerFaceColor', colsubs(M(subjo{ival}),:),'MarkerFaceAlpha', .45,...
            'MarkerEdgeColor', colsubs(M(subjo{ival}),:))

        line([cx-1+rn*.1, cx+rn*.1],...
            [valo(ival), valm(ival)],...
            'Color',...
            [color_lines(dic(sign([valm(ival) - valo(ival)]... determines the color based on the slope of change
            /abs([-(cx-1+rn*.1) + (cx+rn.*1)]))),:), .25],...
            'LineWidth', 1.75)
        hold off
    end
    bc = bc +2;
end
set(ax(1),'XTick', [0, 2],...
    'XTickLabel', {'SE' 'SJ'},...
    'YTick', [-1.5, 0, 4.5], 'LineWidth', 1.5, 'FontName', 'Arial', 'FontSize', 16)

set(ax(2),'XTick', unique(sort([.5:2:3,(.5:2:3)-.5 ,(.5:2:3)+.5])),...
    'XTickLabel', {'OFC', 'SE', 'vmPFC', 'OFC', 'SJ', 'vmPFC'},...
    'YTick', [-.5, 0, 3], 'LineWidth', 1.5, 'FontName', 'Arial', 'FontSize', 16)

ax(1).TickLength(1) = .025;
ax(2).TickLength(1) = .025;
box off % remove the box around the plot
% axes(ax(1)) % select axis 1
sgtitle('\rmvmPFC-OFC -- same brain') % add title to the plot
% print the plot to file
print -dsvg  results\Fig1E.svg
print -dpng -r300  results\Fig1E.png
% write the result to json file
result.p = 2*(1-tcdf(diff(result.mean)./...
    ((diff(result.high) - diff(result.low))/(2*1.96)),length(valm)-1));
result.meanc = diff(result.mean);
result.cic  = result.high-result.low;
result.stdc = (diff(result.high) - diff(result.low))/(2*1.96)*sqrt(length(valm));
json_txt = jsonencode(result, "PrettyPrint",true);
fid = fopen('results\vmPFCOFC_same_brain_stimlock.json', 'w');
fprintf(fid, json_txt);
fclose(fid);


end
% $END
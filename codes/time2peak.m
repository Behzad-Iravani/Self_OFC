%-*- UTF8 -*-
classdef time2peak
    % CLASS NAME: time2peak
    %
    % Purpose: time2peak provides methods for reporting the individual
    % onset latency for HFB self-referential response
    %
    % "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
    %
    % Properties:
    %   - SE: a table containing the summary of data for self-episodic
    %   - SJ: a table containing the summary of data for self-judgment
    % Methods:
    %
    %   Copyright (C)  Behzad Iravani, department of neurology and neurological
    %   sciences, Stanford University. May 2023.
    %
    %   Author: Behzad Iravani
    %   behzadiravani@gmail.com
    %   Contact: behzadiravani@gmail.com
    %   Date: 05/06/2023
    % Author: Behzad Iravani
    % Contact: behzadiravani@gmail.com
    % Date: 05/06/2023
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        SE table % self-episodic data
        SJ table % self-judgment data
    end
    properties(Dependent)
        HFB_response_latencey
    end

    methods
        function obj = time2peak(SE, SJ)
            % time2peak Construct an instance of this class
            obj.SE = SE;
            obj.SJ = SJ;
        end % constructor


        function time2peak = calculate_time2peak(obj)
            %   calculate_time2peak estimates the peak latency for HFB envelope across the
            %   anatomical sites.
            %   Output:
            %      - time2peak: a matrix containing the time (% RT) for the prominent
            %      peak of HFB for every subjects
            %
            %   calculate_time2peak is part of the scripts that recreates the plots
            %   that were reported in "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
            %

            % -----------------------------------------------------------

            Tind      = [obj.SE;obj.SJ]; % concatenates the SE/EP and SJ table
            if iscategorical(Tind.task)% check if the task column is categorical
                Tind.task = cellstr(Tind.task); % converts back to string
            end % if categorical
           s = .1; % smoothing factor

            isubj = 0; % set the subject counter for the for loop
            for subj = unique(Tind.subj)' % loop over all the unique subjects in the data
                isubj = isubj +1; % increment the isubj counter
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
                    % create a time axis as percentage of reaction time
                    time         = linspace(0,100, length(time_series));

                    % find the peak for the smoothed data using matlab findpeaks algorithm
                    [PKs,locs, ~,~] = findpeaks(time_series);
                    % find the peak with the maximum prominence
                    [locs_p]  = find(PKs>quantile(PKs,.75) & time(locs)>10, 1,"first");
                    
                    % storing time ro peak for the give site
                    time2peak.(JP)(isubj) = time(locs(locs_p));
                end % for sites JP
            end % subjs
        end % calculate_time2peak

        function plot(obj)
            dic = containers.Map([-1,0,1],1:3); % a dictionary for converting numbers to index
            color_lines = arrayfun(@(j) misc.hex2rgb(j),["#C60071", "#5E5E5E", "#0073A7"],'UniformOutput',false); % define color for lines   
            % concat the colors
            color_lines = cat(1,color_lines{:});
            % define colors for subjects
            col_subjs = lines(min(...
                            structfun(...
                                @length,(obj.HFB_response_latencey)...
                                )...
                                )); % get colors for subjects
            isubj = size(col_subjs,1); % number of subjects
            hold on % keep the plot 
            rng(10) % for reproducibility
            rn = .15*randn(2,isubj ); % get normally distribute random number of scatter plot
            for isubj = 1:isubj
                % plot OFC data-point for the subject 
                scatter(1+rn(1,isubj), obj.HFB_response_latencey.OFC(isubj), 80,...
                    'filled',...
                    'MarkerFaceColor', col_subjs(isubj,:))
                % plot vmPFC data-point for the subject 
                scatter(2+rn(2,isubj), obj.HFB_response_latencey.MPFC(isubj), 80, ...
                    'filled',...
                    'MarkerFaceColor', col_subjs(isubj,:))
                % plot connecting lines 
                line([1+rn(1,isubj), 2+rn(2,isubj)],...
                    [obj.HFB_response_latencey.OFC(isubj), obj.HFB_response_latencey.MPFC(isubj)],...
                    'Color',...
                    [color_lines(dic(sign([-obj.HFB_response_latencey.OFC(isubj)+ obj.HFB_response_latencey.MPFC(isubj)]... determines the color based on the slope of change
                    /[-1-rn(1,isubj)+2+rn(2,isubj)])),:), .25],...
                    'LineWidth', 1.75)
            end
            % adjusting the plot
            xlim([0,3])
            ylim([0, 100])
            set(gca, 'XTick', 1:2, 'XTickLabel', {'OFC', 'vmPFC'},...
                'FontName', 'Arial Nova Cond', 'FontSize', 24, 'LineWidth', 3,...
                'TickLength', [.02 .01])
            pbaspect([.35,1,1])
            hold off
            print -dsvg results\Fig1c.svg
        end % plot
        % ------------------
        % get methods
        function ROL = get.HFB_response_latencey(obj)
            ROL = obj.calculate_time2peak;
        end
    end % methods
end % class time2peak
% $END
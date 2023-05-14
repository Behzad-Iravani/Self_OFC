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
        mdl(1,1) % mixed-effect model fitted on the group level
    end
    properties(Dependent)
        HFB_response_latencey
    end

    methods
        function obj = time2peak(SE, SJ, mdl)
            % time2peak Construct an instance of this class
            obj.SE = SE;
            obj.SJ = SJ;
            obj.mdl = mdl; 
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
                time2peak.subj(isubj) = subj;
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
                    time = linspace(100*(-.25/max(Tind.time{1}(1,:))),100, length(time_series));

                    % find the peak for the smoothed data using matlab findpeaks algorithm
                     [time2peak.tp.(JP)(isubj), time2peak.p.(JP)(isubj)]= obj.find_peaks(time_series, time, s);   
                     time2peak.rt.(JP)(isubj) = mean(Tind{index,"RT"});
                     time2peak.vrt.(JP)(isubj) = mean(Tind{index,"varRT"});
                     
                end % for sites JP
            end % subjs
          % create within same brain dataset
            ds = table(log10([time2peak.tp.OFC' % log transformation for normality
                time2peak.tp.MPFC']), 0.*[time2peak.p.OFC' 
                time2peak.p.MPFC'], 0.*[time2peak.rt.OFC'
                time2peak.rt.MPFC'], [repmat({'OFC'}, length(time2peak.tp.OFC),1) 
                repmat({'MPFC'}, length(time2peak.tp.MPFC),1)],...
                [repmat(time2peak.subj', 2, 1)],....
                'VariableNames',{'tp', 'p','RT', 'JPAnatomy', 'subj'});
          % control for the variences explined by the model using insample
          % prediction 
          rtp = 10.^(max(ds.tp)).*... rescale to orginal scale
              ((10.^obj.mdl.predict(ds, 'Conditional', true))... invert log10 
              ./max(10.^obj.mdl.predict(ds, 'Conditional', true)));

           time2peak.rtp.OFC = rtp(strcmp(ds.JPAnatomy, 'OFC'));
           time2peak.rtp.MPFC = rtp(strcmp(ds.JPAnatomy, 'MPFC'));
          
        end % calculate_time2peak

        function plot(obj)
            dic = containers.Map([-1,0,1],1:3); % a dictionary for converting numbers to index
            color_lines = arrayfun(@(j) misc.hex2rgb(j),["#C60071", "#5E5E5E", "#0073A7"],'UniformOutput',false); % define color for lines
            % concat the colors
            color_lines = cat(1,color_lines{:});
            % define colors for subjects
            col_subjs = lines(min(...
                structfun(...
                @length,(obj.HFB_response_latencey.rtp)...
                )...
                )); % get colors for subjects
            isubj = size(col_subjs,1); % number of subjects
            hold on % keep the plot
            rng(10) % for reproducibility
            rn = .15*randn(2,isubj ); % get normally distribute random number of scatter plot
            for isubj = 1:isubj
                % plot OFC data-point for the subject
                scatter(1+rn(1,isubj), obj.HFB_response_latencey.rtp.OFC(isubj), 90,...
                    'filled',...
                    'MarkerFaceColor', col_subjs(isubj,:), 'MarkerFaceAlpha', .45,...
                    'MarkerEdgeColor', col_subjs(isubj,:))
                % plot vmPFC data-point for the subject
                scatter(2+rn(2,isubj), obj.HFB_response_latencey.rtp.MPFC(isubj), 90, ...
                    'filled',...
                    'MarkerFaceColor', col_subjs(isubj,:), 'MarkerFaceAlpha', .45,...
                    'MarkerEdgeColor', col_subjs(isubj,:))
                % plot connecting lines
                line([1+rn(1,isubj), 2+rn(2,isubj)],...
                    [obj.HFB_response_latencey.rtp.OFC(isubj), obj.HFB_response_latencey.rtp.MPFC(isubj)],...
                    'Color',...
                    [color_lines(dic(sign([-obj.HFB_response_latencey.rtp.OFC(isubj)+ obj.HFB_response_latencey.rtp.MPFC(isubj)]... determines the color based on the slope of change
                    /[-1-rn(1,isubj)+2+rn(2,isubj)])),:), .25],...
                    'LineWidth', 1.75)
            end
            % adjusting the plot
            xlim([0,3])
            ylim([0, 100])
            set(gca, 'XTick', 1:2, 'XTickLabel', {'OFC', 'vmPFC'},...
                'FontName', 'Arial Nova Cond', 'FontSize', 24, 'LineWidth', 3,...
                'TickLength', [.02 .01])
            ylabel('RT%')
            title({'\rm LMM in-sample prediction', 'individuals with OFC & vmPFC'}, 'FontName', 'Arial Nova Cond', 'FontSize', 12)
            pbaspect([.35,1,1])
            hold off
            print -dsvg results\Fig1c_2.svg
        end % plot
        % ------------------
        % get methods
        function ROL = get.HFB_response_latencey(obj)
            ROL = obj.calculate_time2peak;
        end
    end % methods
    methods(Static)
        function [t, p]= find_peaks(S, T, s)
            S = movmean(S, floor(s.*1e3));
            if sum(T<0)
            S = S -mean(S(T<0));
            end
            % find the peak for the smoothed data using matlab findpeaks algorithm
            [PKs,locs, ~,~] = findpeaks(S, 'MinPeakWidth', .1, 'MinPeakDistance', .1);
            % find the peak with the maximum prominence
            [locs_p]  = find(PKs>0 & PKs>quantile(PKs,.975) & T(locs)>10, 1,"first");
            
            t = T(locs(locs_p));
            if isempty(t)
                t = nan;
            end
            p = PKs(locs_p);
            if isempty(p)
                p = nan;
            end
        end % find_peak

    end % static
end % class time2peak
% $END
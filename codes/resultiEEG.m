%-*- UTF8 -*-
classdef resultiEEG < stat_report
    % CLASS NAME: resutiEEG
    %
    % Purpose: resutiEEG provides methods for recreating the results and figures reported in
    % "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX ",
    %
    %
    % Properties:
    %  subclass of stat_report
    %
    %
    % Methods:
    %     LocalizeSelfMath: plots the self-referentially vs. math activated electrodes on the MNI brain surface
    %
    % Author: Behzad Iravani
    % Contact: behzadiravani@gmail.com
    % Date: 05/02/2023
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    methods
        function obj = resultiEEG(data, jbpath, jepath)
            % constructor method for creating instance of resultiEEG
            obj@stat_report(data, jbpath, jepath)
        end

        function obj = LocalizeSelfMath(obj, col)
            % LocalizeSelfMath plots electrodes activity on the brain surface
            % Input:
            %      - col:   colors used for mapping electrode activities
            %-----------------------------------------------------------
            % Regroup conditions (EP: self-episodic and SJ: self-judgment)
            % as self.
            sub = obj.localizing(obj.data);
           
            % divide to left and right self-referential data
            right.self    = select_hemi(sub{1}, "right");
            left.self     = select_hemi(sub{1}, "left");

            % adding elec type
            right.self    = adding_electype(right.self);
            left.self     = adding_electype(left.self);

            % divide to left and right math data
            right.mth     = select_hemi(sub{2}, "right");
            left.mth      = select_hemi(sub{2}, "left");

            % adding elec type
            right.mth     = adding_electype(right.mth);
            left.mth      = adding_electype(left.mth);


            c = 0; % sets the counter for plotting
            for hemi = {'right', 'left'} % loop over the hemispheres
                clear MNI_coord DATA
                c = c+1; % increment the counter
                if strcmp(hemi{:}, 'right')
                    MNI_coord = [right.self.X(right.self.X>0), right.self.Y(right.self.X>0), right.self.Z(right.self.X>0)];
                    DATA = right.self.Loc_Tval; % localizing stats Self > Math
                    % remove those that were de-activated during math
                    DATA(right.mth.Tval<0 & right.mth.Pval<0.05) = 0;
                else % left
                    MNI_coord = [left.self.X(left.self.X<0), left.self.Y(left.self.X<0), left.self.Z(left.self.X<0)];
                    DATA = left.self.Loc_Tval; % localizing stats Self > Math
                    % remove those that were de-activated during math
                    DATA(left.mth.Tval<0 & left.mth.Pval<0.05) = 0;
                end
                if ~any(DATA>=0) % checks if data exists
                    continue
                end
                % create surface
                surface = obj.Create_surface(MNI_coord, DATA);
                % plot coverage
                obj.plot_coverage(surface, col, hemi, right.self, left.self)
                % add normalized left/right T-values
                right.self.Tn = right.self.Loc_Tval/max(abs([right.self.Loc_Tval; left.self.Loc_Tval]));
                left.self.Tn  = left.self.Loc_Tval/max(abs([right.self.Loc_Tval; left.self.Loc_Tval]));
                % plot electrodes on the brain
                plot_electActivity(surface, hemi, right.self, left.self, col)
            end % for

            function subT = select_hemi(T, hemi)
                % select_hemi divides the table to two tables containing the right and left hemispheres data
                % Input:
                %      T: table containing data
                %      hemi: string defines left or right
                % ---------------------------------
                switch hemi % determines if left or right
                    case "right"
                        subT = T(T.X>0, :);
                    case "left"
                        subT = T(T.X<0, :);
                    otherwise
                        error("hemi should be either left or right.")
                end
                ind(:,1) = subT.Loc_Tval>0 & subT.Loc_Pval <= .05; %  more than math
                ind(:,2) = subT.Loc_Tval<0 & subT.Loc_Pval <= .05; % more than self

                subT.select(ind(:,1)) = {'Self'};
                subT.select(ind(:,2)) = {'Math'};
                subT.select(~ind(:,1) & ~ind(:,2) ) = {'none'};
                subT.select = categorical(subT.select);
            end % select_hemi

            function Tout = adding_electype(Tin)
                Tout = Tin;
                for s = unique(Tout.subj)'
                    [~, ~, ib] = intersect(s, {obj.ECoGSEEG.subj});
                    Tout.elecType(strcmp(Tout.subj, s{:})) =  {obj.ECoGSEEG(ib).elec_type}';
                end % for subj
            end % adding_electype
        end % LocalizeSelfMath

    end % methods
    methods(Static)
        function surface = Create_surface(MNI, data)
            %  Static method for surface object
            % Input
            surface  = surf_();

            surface.left.path      ='codes\@surf_\dat\lh.pial';
            surface.right.path     ='codes\@surf_\dat\rh.pial';
            surface.ElectPos       = MNI;
            surface.ElectActivity  = data;

        end

        function HFB = getTimeWarppedHFB(file_)
            % Reads time warped HFB enveloped by the EEGLAB time-warp
            % function from *.dat file
            % Input:
            %       file_:  string or character vector containing the path to *.dat file
            %
            % Output:
            %       HFB: matrix of time-warped HFB, labels exists in corresponding *.tsv file
            %---------------------------------------------------------

            % Open the .dat file
            fid = fopen(file_, 'r');

            % Read the size of the matrix from the file
            matrix_size = fread(fid, [1 2], 'int');

            % Read the matrix from the file
            HFB = fread(fid, fliplr(matrix_size), 'double')';

            % Close the file
            fclose(fid);

        end % getTimeWarpedHFB

        function mdl1 = plot_HFB(HFB, s, col)
            % plot_HFB plots the HFB time course for the anatomical sites
            % Input:
            %       - HFB: structure with filed:
            %              - data: matrix containing the
            %                time warped HFB enveloped
            %               - label:
            %                 table containing the labels for data
            %       - s:  scalar, smoothing factor for visualization purposes.
            %       - col:
            %             string array containing the hex color codes.
            %-------------------------------------------------------------
            iplot = 0;
            figure
            clear p
            hold on
            fsamp = round(1/median(diff(HFB.time(1,:))));
            clear T
            for JP = ["MPFC","OFC"]
                iplot =  iplot +1;
                %     subplot(1,2,iplot)
                ci = 0; % set the counter
                for tsk = ["self"] % for self-referential task
                    index = cellfun(@(x) any(strcmp(x, JP)), HFB.label.JPAnatomy); % find the indices for the given sites 
                    ci = ci + 1; % increment the counter
                    % assemble the table for the given site and
                    % self-referential task
                    T{iplot} = HFB.label(index,:);
                    % select 250 to end
                    T{iplot}.time = num2cell(HFB.time(index, sum(HFB.time>-.250)>0),2); % add time to the table
                    % add the ieeg data (HFB)
                    T{iplot}.avg = num2cell(HFB.data(index,sum(HFB.time>-.250)>0),2);
                    % average over the task
                    T{iplot} =  resultiEEG.localizing(T{iplot});
                    T{iplot}{1} = T{iplot}{1}(T{iplot}{1}.self_actiavted == 1,:);
                    [time_series, all.(JP), nn.(JP), ss.(JP), ci_.(JP)] = misc.t_varaible_length(T{iplot}{1},'stim', s); % computes t-value for every time bins
                    tt.(JP) = movmean(time_series,floor(s*1e3));
                    time = linspace(100*(-.25/max(T{iplot}{1}.time{1})),100, length(tt.(JP))); % converting the time to percentage (0 to 100%) of RT
                    % plot HFB traces
                    p(iplot) = plot(time, movmean(tt.(JP),floor(s*1e3)), 'color', misc.hex2rgb(col(iplot)), 'LineWidth', 2);

                    % add shadings for confidence intervals
                    fill([time, fliplr(time)], [movmean(tt.(JP),floor(s*1e3))-1.96, fliplr(movmean(tt.(JP),floor(s*1e3)))+1.96], misc.hex2rgb(col(iplot)),...
                        'FaceAlpha', .12, 'EdgeColor', 'none' )

                    ylim([-15,15])
                    set(gca, 'FontName', 'Arial', 'FontSize', 20, 'LineWidth', 1.75, 'YTick', [-15:15:15])
                    line(xlim(), [2 2], 'LineStyle', '--', 'Color', 'k')
                    line(xlim(), [-2 -2], 'LineStyle', '--', 'Color', 'k')
                end % for tsk
                box off
                xlabel('RT%')
                ylabel('HFB (t-value)')
            end % for sites -- JPAnatomy
            axis square % make the axis square


            % write the result in the *.json file
            [~, index_diff] = findpeaks(double(tt.MPFC-tt.OFC>(tinv(.99, nn.MPFC) + tinv(.99,nn.OFC))/2 & ...
                tt.MPFC>0 & tt.OFC>0),...
                'MinPeakDistance', round(.1*fsamp),... fsamp = 512
                'MinPeakWidth',round(.1*fsamp)...
                ); % find the index that vmPFC is larger than OFC

            [~, index_MPFC] = findpeaks(double(tt.MPFC>(tinv(.99, nn.MPFC))),...
                'MinPeakDistance', round(.1*fsamp),... fsamp = 512
                'MinPeakWidth',round(.1*fsamp)...
                ); % find the index that vmPFC is larger the critical t
            [~, index_OFC] = findpeaks(double(tt.OFC>(tinv(.99, nn.OFC))),...
                'MinPeakDistance', round(.1*fsamp),... fsamp = 512
                'MinPeakWidth',round(.1*fsamp)...
                ); % find the index that OFC is larger the critical t
           % run a control LMM analysis
           [T{1}{1}.tp, T{1}{1}.p] = cellfun(@(x)time2peak.find_peaks(x, time,.1), T{1}{1}.avg);
           [T{2}{1}.tp, T{2}{1}.p] = cellfun(@(x)time2peak.find_peaks(x, time,.1), T{2}{1}.avg);
           % concatenate tables and rum LMM
            TL = [T{1}{1}; T{2}{1}];
            % log transformation for normality 
            TL.RT = log10(TL.RT);
            TL.tp = log10(TL.tp);
            
           
            mdl1 = fitlme(TL, ...
                ['tp ~ 1 + JPAnatomy * p + RT +'....
                '(1 + JPAnatomy |subj) '], 'DummyVarCoding','reference', ...
                'Verbose',false);

            A = anova(mdl1);
            result.mdl.Estimate = mdl1.Coefficients.Estimate;
            result.mdl.SE     = mdl1.Coefficients.SE;
            result.mdl.DF     = mdl1.Coefficients.DF;
            result.mdl.tStat  = mdl1.Coefficients.tStat;
            result.mdl.pValue = mdl1.Coefficients.pValue;
            result.mdl.Lower  = mdl1.Coefficients.Lower;
            result.mdl.Upper = mdl1.Coefficients.Upper;
                
            result.anova.term = A.Term;
            result.anova.FStat = A.FStat;
            result.mdl.DF1 = A.DF1;
            result.mdl.DF2 = A.DF2;
            result.mdl.pValue = A.pValue;
           
            result.Anatomy = fieldnames(tt);
            result.timeRelative_percentage    = round([time(index_diff(1))]);
            result.timeAbsolute_percentage    = round([time(index_MPFC(1)), time(index_OFC(1))]);
            result.tvalue  = [tt.MPFC(index_diff(1)), tt.OFC(index_diff(1))];
            result.dof     = [nn.MPFC(index_diff(1))-1, nn.OFC(index_diff(1))-1];
            result.pvalue  = 2*(1-tcdf(abs(result.tvalue), result.dof));
            result.CI      = [ci_.MPFC(:,index_diff(1))'
                ci_.OFC(:,index_diff(1))'];

            json_text = jsonencode(result, "PrettyPrint",true); % encode the struct data to json text
            disp(json_text) % print the text to console
            fid = fopen('results\vmPFCOFC_faster.json', 'w'); % create and open the file
            fprintf(fid, json_text); % writes the text to file
            fclose(fid); % close the file
            scatter(time(index_diff(1)), 1.10*tt.MPFC(index_diff(1)), 'filled', 'v', 'MarkerFaceColor', 'r', 'MarkerEdgeColor','none')
            % add legends
            legend(p, {'vmPFC', 'OFC'}, 'Box','off', 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontName', 'Arial')

            hold off
            % printing results
            print results\Fig1c.svg -dsvg
            print results\Fig1c.png -r300 -dpng

        end %  plot_HFB
       

        function plotHrdWr(HFB, s, col)

            % SEEG
            index.SEEG = (strcmp(HFB.label.subj, 'S17') | ...
                strcmp(HFB.label.subj, 'S19') | ...
                strcmp(HFB.label.subj, 'S21') | ...
                strcmp(HFB.label.subj, 'S22') )& ...
                (strcmp(HFB.label.task, 'EP') | ...
                strcmp(HFB.label.task, 'SJ') )& ...
                (HFB.label.Loc_Tval>2 & HFB.label.Loc_Pval<0.05 & HFB.label.Tval>.5);

            index.ECoG =~(strcmp(HFB.label.subj, 'S17') | ...
                strcmp(HFB.label.subj, 'S19') | ...
                strcmp(HFB.label.subj, 'S21') | ...
                strcmp(HFB.label.subj, 'S22') )& ...
                (strcmp(HFB.label.task, 'EP') | ...
                strcmp(HFB.label.task, 'SJ') )& ...
                (HFB.label.Loc_Tval>2 & HFB.label.Loc_Pval<0.05 & HFB.label.Tval>.5);
            % find the closest SEEG and ECoG
            SEEGECOGDIS = sqrt((HFB.label.X(index.SEEG) - HFB.label.X(index.ECoG)').^2 +...
                (HFB.label.Y(index.SEEG) - HFB.label.Y(index.ECoG)').^2 +...
                (HFB.label.Z(index.SEEG) - HFB.label.Z(index.ECoG)').^2 );
            [~, index.c] = min(SEEGECOGDIS,[],'all');
            [index.cS(1), index.cS(2)] = ind2sub(size(SEEGECOGDIS), index.c);
            index.SEEG = find(index.SEEG);
            index.ECoG = find(index.ECoG);

           index.SEEG = find( strcmp(HFB.label.subj,HFB.label.subj(index.SEEG(index.cS(1)))) &...
                        strcmp(HFB.label.chan, HFB.label.chan(index.SEEG(index.cS(1)))) & ...
                        (strcmp(HFB.label.task, 'EP') | ...
                        strcmp(HFB.label.task, 'SJ')));

           index.ECoG = find( strcmp(HFB.label.subj,HFB.label.subj(index.ECoG(index.cS(2)))) &...
                        strcmp(HFB.label.chan, HFB.label.chan(index.ECoG(index.cS(2)))) & ...
                        (strcmp(HFB.label.task, 'EP') | ...
                        strcmp(HFB.label.task, 'SJ')));

            iplot = 0;
            % converting the time to percentage (0 to 100%) of RT
            hold on
            for p = ["SEEG", "ECoG"]
                iplot = iplot +1;
                T = HFB.label(index.(p),:);
                % select 250 to end
                T.time = num2cell(HFB.time(index.(p), sum(HFB.time>-.250)>0),2); % add time to the table
                % add the ieeg data (HFB)
                T.avg = num2cell(HFB.data(index.(p),sum(HFB.time>-.250)>0),2);

                [time_series, all.(p), nn.(p), ss.(p), ci_.(p)] = misc.t_varaible_length(T,'stim', s); % computes t-value for every time bins
                time = linspace(100*(-.25/max(T.time{1})),100, size(all.(p),2));
                tt.(p).mean = mean(movmean(time_series(T.time{1}>.8 & T.time{1}<1.8),floor(s*1e3)));
              
                % plot HFB traces

                tt.(p).mtval = max(T.Tval); %tt.(p)(T.time{1} > .8 & T.time{1} < 1.8);
                tt.(p).dof = fix(mean(T.dof));
                tt.(p).SE  =  (tt.(p).mean./tt.(p).mtval);
                tt.(p).var = (tt.(p).dof+1).*(tt.(p).SE).^2; 
                pl{iplot} = bar(iplot, tt.(p).mean ,...
                    'FaceColor', misc.hex2rgb(col(iplot))');
                errorbar(iplot, tt.(p).mean, tt.(p).SE ...std(Sig)
                    ,'Color', 'k', 'LineWidth', 1.75, 'CapSize', 0)
                %                     plot(time, zscore(),...
                %                     'color', misc.hex2rgb(col(iplot)), 'LineWidth', 2);

                % add shadings for confidence intervals
                %                 fill([time, fliplr(time)], [zscore(movmean(time_series.(p),floor(s*1e3)))-1.96,...
                %                     fliplr(zscore(movmean(time_series.(p),floor(s*1e3)))+1.96)],...
                %                     misc.hex2rgb(col(iplot)),...
                %                     'FaceAlpha', .15, 'EdgeColor', 'none' )


                %                 for iplot2 = 1:size(all.(p))
                %                     pl{iplot}(iplot2) = plot(time, movmean(zscore(all.(p)(iplot2,:)),floor(s*1e3)), 'color', misc.hex2rgb(col(iplot)), 'LineWidth', 2);
                %
                %                     % add shadings for confidence intervals
                %                     fill([time, fliplr(time)], [movmean(zscore(all.(p)(iplot2,:)),floor(s*1e3))-1.96,...
                %                         fliplr(movmean(zscore(all.(p)(iplot2,:)),floor(s*1e3))+1.96)],...
                %                         misc.hex2rgb(col(iplot)),...
                %                         'FaceAlpha', .05, 'EdgeColor', 'none' )
                %                 end % for iplot2
            end % for p
            %             uistack(pl{2}, 'top'); % bring the lines to top
            %             uistack(pl{1}, 'top')
            sprintf('t(%d) = %1.2f , p = %1.2f\n', (tt.SEEG.dof + tt.ECoG.dof), abs(tt.SEEG.mean - tt.ECoG.mean)./sqrt(tt.SEEG.var + tt.ECoG.var),...
                2*(1-tcdf(abs(tt.SEEG.mean - tt.ECoG.mean)./sqrt(tt.SEEG.var + tt.ECoG.var), (tt.SEEG.dof + tt.ECoG.dof))))
            box off
            set(gca, 'XTick', [1,2], ...
                'XTickLabel', {'SEEG', 'ECoG'},'FontName', 'Arial Nova Cond', 'FontSize', 18, 'LineWidth', 2)

            ylabel('Normalized HFA power (STD)')
            %             legend([pl{1}, pl{2}], ...
            %                 arrayfun(@(x)sprintf('n = %d',x ),[median(nn.SEEG), median(nn.ECoG)],'UniformOutput',false),...
            %                 'FontName', 'Arial Nova Cond', 'Box', 'off')
           xlim([0,3])
           ylim([0,.05])
            pbaspect([.5,1,1])
            print -dpng -r300 results\figs2b_2.png
            print -dsvg -vector results\figs2b_2.svg

         
            
%              resultiEEG.muI(tt.SEEG(T.time{1} > .8 & T.time{1} < 1.8)', tt.ECoG(T.time{1} > .8 & T.time{1} < 1.8)', ...
%                 {linspace(min([tt.SEEG(T.time{1} > .8 & T.time{1} < 1.8)'
%                 tt.ECoG(T.time{1} > .8 & T.time{1} < 1.8)']), ...
%                 max([tt.SEEG(T.time{1} > .8 & T.time{1} < 1.8)' 
%                 tt.ECoG(T.time{1} > .8 & T.time{1} < 1.8)']), 10),...
%                 linspace(min([tt.SEEG(T.time{1} > .8 & T.time{1} < 1.8)' 
%                 tt.ECoG(T.time{1} > .8 & T.time{1} < 1.8)']), ...
%                 max([tt.SEEG(T.time{1} > .8 & T.time{1} < 1.8)'
%                 tt.ECoG(T.time{1} > .8 & T.time{1} < 1.8)']), 10)}...
%                 )
        end % plotHrdWr

        function [EP, SJ] = find_sameBrain(T, tcs, time)
            % find_sameBrain finds the patients with electrodes that were implanted in
            % both OFC and vmPFC
            % Input:
            %       -T:     a table containing the data regrading subjects,
            %       channels and etc.
            %       -tcs:   a matrix containing HFB time-course data
            %       -time:  a matrix containing time axis
            %--------------------------------------------------
            if nargin>1
                T.avg = num2cell(tcs,2); % converting the matrix to cell array.
                T.time = repmat({time}, height(T),1); % repate the time axis to a cell array.
            end
            if iscategorical(T.subj) % checks if the subj column is categorical
                T.subj =  cellstr(T.subj); % convert to string if it is categorical
            end % if iscategorical subj
            if ~iscategorical(T.subj) % checks if the task column is not categorical
                T.task =  categorical(T.task); % convert to categorical if it is string
            end % if iscatehorical task
            subs = []; % initialize the subj variable that stores the subject name
            EP   = []; % initialize the EP variable that stores the EP(self-episodic/SE) data
            SJ   = []; % initialize the SJ variable that stores the SJ(self-judgment) data

            for s = unique(T.subj)' % loop over the unique subject names
                subs = [subs;s]; % store the new subject's name in variable subj
                if any(T.task == 'EP') % checks if the task is self-episodic
                    l = numel(unique(T.JPAnatomy(strcmp(T.subj, s{:})& ...
                        T.task == 'EP'))); % find out if more than one datapoint exists
                    if l>1 % checks if more than one datapoint exists -> average over and store in EP
                        EP = [EP; stat.average_over_rois(T(strcmp(T.subj, s{:}) & ...
                            T.task == 'EP',:))];
                    end% if l
                end % if EP(SE)
                if any(T.task == 'SJ') % checks if the task is self-judgment
                    l = numel(unique(T.JPAnatomy(strcmp(T.subj, s{:})& ...
                        T.task == 'SJ'))); % find out if more than one datapoint exists
                    if l>1  % checks if more than one datapoint exists -> average over and store in SJ
                        SJ = [SJ; stat.average_over_rois(T(strcmp(T.subj, s{:}) & ...
                            T.task == 'SJ',:))];
                    end% if l
                end % if SJ
            end % for
        end % sameBrain
        function plot_coverage(surface, col, hemi, right, left)
            % plot coverage
            % load MNI surface of hemi
            [vertex_coords, faces] = surface.read_surf(surface.(hemi{:}).path);

            if ~isempty(find(faces == 0)) % checks if the indexing starts with 0
                warning('indexing starts at zero adding 1 to faces')
                faces = faces +1;
            end
            % assessing the vertexes and face to surface object
            surface.(hemi{:}).cortex.vert  = vertex_coords;
            surface.(hemi{:}).cortex.tri   = faces;
            % computes density map and correcting
            MNI_surface = surface.densityMAP(hemi{:}); % computes the density of electrodes
            for views = {'ventral', 'medial'}
                figure % open a figure
                % plot brain surface
                [ax, trs]  = surface.plot_brain(hemi, [.95 .95 .85]);
                alpha(1)
                hold on
                if strcmp(hemi{:}, 'right')
                    actvie_elec = right.select ;
                    p = right.Tval;
                    [s, D, flag, label] = misc.plot_electrodes_projected_to_brain(...
                        [right.X, right.Y, right.Z], surface, hemi, views{:}, right.JPAnatomy, 25);
                else
                    actvie_elec = left.select ;
                    p = left.Tval;
                    [s, D, flag, label] = misc.plot_electrodes_projected_to_brain(...
                        [left.X, left.Y, left.Z], MNI_surface, hemi, views{:}, left.JPAnatomy, 25);
                end
                misc.reorgnize_electrodes(s, D, actvie_elec, flag, col, right, left, hemi, views, true)
                misc.plot_setting(hemi{:},views{:})
                axis equal
                axis tight
                axis off
                camlight
                print('-dpng', '-r300', ['results\figS2b_', hemi{:} '_' views{:}, '.png'])
            end % views
        end % plot coverage
        function mutualInfo = muI(x,y, bins)
            jointHist = hist3([x y], 'Edges', bins);
            % Normalize the joint histogram to obtain the joint probability distribution
            jointProb = jointHist / numel(x);
            figure
            contourf(bins{1}, bins{2}, jointProb,30, 'LineStyle', 'none')
            ylabel('SEEG')
            xlabel('ECoG')
            axis square
            set(gca, ...
                'FontName', 'Arial Nova Cond', 'FontSize', 18, 'LineWidth', 2)
            colormap(flipud(gray(32)))
            title('\rm Joint distribution')
            box off
            % Compute the marginal probability distributions
            marginalX = sum(jointProb, 2);
            marginalY = sum(jointProb, 1);

            % Compute the entropy of x
            entropyX = -nansum(marginalX .* log2(marginalX));

            % Compute the entropy of y
            entropyY = -nansum(marginalY .* log2(marginalY));

            % Compute the joint entropy
            entropyJoint = -nansum(jointProb(:) .* log2(jointProb(:)));

            % Compute the mutual information
            mutualInfo = entropyX + entropyY - entropyJoint;

            text(.5,.5,['Mutual Information: ' num2str(round(mutualInfo),2)],...
                'Unit', 'Normalized', 'Color', 'k', ...
                'FontName', 'Arial Nova Cond', 'FontSize', 12,'HorizontalAlignment', 'Center');
            print -dsvg -vector results\figs1b3.svg
        end
        function sub = localizing(Tin)
            % localizing find self-referntially active electrodes that
            % were active during self-referential but not de-active in math
            % average over self 
              % Regroup conditions (EP: self-episodic and SJ: self-judgment)
            % as self.
            itsk = 0; % task counter
            for tsk = ["self", "MTH"] % first and second elements of sub are self and math respectively
                itsk = itsk +1;
                if strcmp(tsk, "MTH")
                    index = strcmp(Tin.task, tsk);
                else
                    index = (strcmp(Tin.task, "EP") | strcmp(Tin.task, "SJ")); % finds both Self-episodic (EP) and self-judgment
                end
                sub{itsk} = Tin(index,:);
                sub{itsk} = stat.average_over_sessions(sub{itsk}); % average over sessions.
            end
            % average over session for the first element that contains two
            % conditions
            sub{1} = stat.average_over_task(sub{1});% first average over task (EP & SJ)
            sub{1}.task = repmat({'self'}, height(sub{1}), 1); % reassign the task to "self"
            % math
            sub{2} = stat.average_over_task(sub{2});% average over MTH

            % find self activated
            sub{1}.self_actiavted = sub{1}.Loc_Tval > 0 & sub{1}.Loc_Pval < 0.05 ... more in self than math
                & ~(sub{2}.Tval < 0 & sub{2}.Pval < 0.05); % not deactivated in math  
            if sum(sub{1}.self_actiavted) == 0
                error('no self-activated electrode exists!')
            end
 
        end

    end % methods static
end % class resultiEEG
% $END
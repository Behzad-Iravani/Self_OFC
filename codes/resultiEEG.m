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
            itsk = 0; % task counter
            for tsk = ["self", "MTH"] % first and second elements of sub are self and math respectively
                itsk = itsk +1;
                if strcmp(tsk, "MTH")
                    index = strcmp(obj.data.task, tsk);
                else
                    index = (strcmp(obj.data.task, "EP") | strcmp(obj.data.task, "SJ")); % finds both Self-episodic (EP) and self-judgment
                end
                sub{itsk} = obj.data(index,:);
                sub{itsk} = stat.average_over_sessions(sub{itsk}); % average over sessions.
            end
            % average over session for the first element that contains two
            % conditions
            sub{1} = stat.average_over_task(sub{1});% first averge over task (EP & SJ)
            sub{1}.task = repmat({'self'}, height(sub{1}), 1); % reassign the task to "self"
            % divide to left and right self-referential data
            right    = select_hemi(sub{1}, "right");
            left     = select_hemi(sub{1}, "left");

            figure % open a new figure
            c = 0; % sets the counter for plotting
            for hemi = {'right', 'left'} % loop over the hemispheres
                clear MNI_coord DATA
                c = c+1; % increament the counter
                if strcmp(hemi{:}, 'right')
                    MNI_coord = [right.X(right.X>0), right.Y(right.X>0), right.Z(right.X>0)];
                    DATA = right.Loc_Tval; % localizing stats Self > Math
                else % left
                    MNI_coord = [left.X(left.X<0), left.Y(left.X<0), left.Z(left.X<0)];
                    DATA = left.Loc_Tval; % localizing stats Self > Math
                end
                if ~any(DATA>=0) % checks if data exists
                    continue
                end
                % create surface
                surface = obj.Create_surface(MNI_coord, DATA);

                % add normalized left/right T-values
                right.Tn = right.Loc_Tval/max(abs([right.Loc_Tval; left.Loc_Tval]));
                left.Tn  = left.Loc_Tval/max(abs([right.Loc_Tval; left.Loc_Tval]));
                % plot electrodes on the brain
                plot_electActivity(surface, hemi, right, left, col)
            end

            function subT = select_hemi(T, hemi)
                % select_hemi divides the table to two tables containing the right and left hemispheres data
                % Input:
                %      T: table containing data
                %      hemi: string defines left or right
                % ---------------------------------
                switch hemi % determines if lef or right
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
            % Reads timewarped HFB enveloped by the EEGLAB timewarp
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
            HFB = fread(fid, matrix_size, 'double');

            % Close the file
            fclose(fid);

        end % getTimeWarpedHFB

        function plot_HFB(HFB, s, col)
            % plot_HFB plots the HFB time course for the anatomical sites
            % Input:
            %       - HFB: structure with filed:
            %              - data: matrix containing the
            %                time warped HFB enveloped
            %               - label:
            %                 table containing the labels for data
            %       - s:  scalar, smoothing factor for visulaiztion purposes.
            %       - col:
            %             string array containing the hex color codes.
            %-------------------------------------------------------------
            iplot = 0;
            figure
            clear p
            hold on
            for JP = ["MPFC","OFC"]
                iplot =  iplot +1;
                %     subplot(1,2,iplot)
                ci = 0; % set the counter
                for tsk = ["self"] % for self-referential task
                    index = cellfun(@(x) any(strcmp(x, JP)), HFB.label.JPAnatomy) & (strcmp(HFB.label.task, "EP") | strcmp(HFB.label.task, "SJ")); % find the indices for the given sites and self-referental tasks (self-episodic: EP, self-judgmen: SJ)
                    ci = ci + 1; % increament the counter
                    % assembel the table for the given site and
                    % self-referntial task
                    T = HFB.label(index,:);
                    T.time = repmat({HFB.time}, height(T), 1); % add time to the table
                    % add the ieeg data (HFB)
                    T.avg = num2cell(HFB.data(index,:),2);
                    % average over the task
                    [time_series, all.(JP), nn.(JP), ss.(JP), ci_.(JP)] = misc.t_varaible_length(T,'stim', s); % computes t-value for every time bins
                    tt.(JP) = movmean(time_series,floor(s*1e3));
                    time = linspace(0,100, length(HFB.time)); % converting the time to precentage (0 to 100%) of RT
                    % plot HFB traces
                    p(iplot) = plot(time, movmean(time_series,floor(s*1e3)), 'color', misc.hex2rgb(col(iplot)), 'LineWidth', 2);

                    % add shadings for confidence intervals
                    fill([time, fliplr(time)], [movmean(time_series,floor(s*1e3))-1.96, fliplr(movmean(time_series,floor(s*1e3)))+1.96], misc.hex2rgb(col(iplot)),...
                        'FaceAlpha', .12, 'EdgeColor', 'none' )

                    ylim([-6,6])
                    set(gca, 'FontName', 'Arial', 'FontSize', 20, 'LineWidth', 1.75, 'YTick', [-6:6:6])
                    line(xlim(), [2 2], 'LineStyle', '--', 'Color', 'k')
                    line(xlim(), [-2 -2], 'LineStyle', '--', 'Color', 'k')
                end % for tsk
                box off
                xlabel('%RT')
                ylabel('HFB (t-value)')
            end % for sites -- JPAnatomy
            axis square

            legend(p, {'vmPFC', 'OFC'}, 'Box','off', 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontName', 'Arial')
            hold off
            % printing results
            print results\Fig1c.svg -dsvg
            print results\Fig1c.png -r300 -dpng

            result.Anatomy = fieldnames(tt);
            result.time    = round([time(find(tt.MPFC>2.05, 1, "first")), time(find(tt.OFC>2.05, 1, "first"))]);
            result.tvalue  = [tt.MPFC(find(tt.MPFC>2.05, 1, "first"))', tt.OFC(find(tt.OFC>2.05, 1, "first"))];
            result.dof     = [nn.MPFC(find(tt.MPFC>2.05, 1, "first"))'-1, nn.OFC(find(tt.OFC>2.05, 1, "first"))-1];
            result.pvalue  = 2*(1-tcdf(abs(result.tvalue), result.dof));
            result.CI      = [ci_.MPFC(:,find(tt.MPFC>2.05, 1, "first"))'
                ci_.OFC(:,find(tt.OFC>2.05, 1, "first"))'];

            json_text = jsonencode(result, "PrettyPrint",true);
            disp(json_text)
            fid = fopen('results\vmPFCOFC_faster.json', 'w');
            fprintf(fid, json_text);
            fclose(fid)

        end %  plot_HFB


    end
end

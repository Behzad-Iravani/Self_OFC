classdef LMMPoNe < stat.LMM
    %   CLASS NAME: LMMSPoNe
    %
    %   Purpose: Linear mixed-effect model for comparing positive vs. negative connotations in
    %   OFC and vmPFC
    %   This script is part of the codes for reproducing stats reported in
    %   "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
    %
    %   Properties:
    %       - index:        a structure containing the indices for each coefficient
    %   Methods:
    %
    %   Copyright (C)  Behzad Iravani, department of neurology and neurological
    %   sciences, Stanford University. May 2023.
    %
    %   Author: Behzad Iravani
    %   behzadiravani@gmail.com
    %   Contact: behzadiravani@gmail.com
    %   Date: 05/09/2023
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        index  % a structure containing the indices for each coefficient
    end
    properties(Dependent)
        preprocT % preprocessed data for True(self-coherent) and False (self-incoherent) LMM
    end
    methods
        function obj = violin(obj)
            vdat = stat.average_over_sessions(obj.preprocT);
            Po =  stat.average_over_task(vdat(strcmp(vdat.task,'one'),:)); % Positive data
            Ne =  stat.average_over_task(vdat(strcmp(vdat.task,'minusone'),:)); % Negative data
            dic = containers.Map([-1,0,1],1:3); % a dictionary that maps the slope to colors index

            col = [0.7753         0    0.4432
                0.3695    0.3695    0.3695
                0    0.4508    0.6564];
            col_pos_neg = obj.colors([1,2],1); % get colors for violin plots
            i = 0;
            for Anat = {'OFC', 'MPFC'}

                % find unique data points for Po and Ne
                [ip, in] = find(categorical(cellfun(@(x,y) [x ':' y],...
                    Po{Po.JPAnatomy == Anat{:},"subj"},...
                    Po{Po.JPAnatomy == Anat{:},"chan"}, 'UniformOutput', false)) ...
                    == ...
                    categorical(cellfun(@(x,y) [x ':' y],...
                    Ne{Ne.JPAnatomy == Anat{:},"subj"},...
                    Ne{Ne.JPAnatomy == Anat{:},"chan"}, 'UniformOutput', false))');

                clear v x_tmp y_tmp
                hold on
                i = i + 1;
                v = Violin(Po.Tval(ip), i);
                % adjust the violin plot
                misc.violin_setting(v)
                % change colors
                v.ScatterPlot.MarkerFaceColor = col_pos_neg(1,:);
                v.ScatterPlot.MarkerEdgeColor = col_pos_neg(1,:);
                v.ViolinPlot.FaceColor        = col_pos_neg(1,:);
                % save coordinates for plotting lines
                x_tmp(1,:) = v.ScatterPlot.XData;
                y_tmp(1,:) = v.ScatterPlot.YData;
                % --------------------------------
                i = i + 1;
                v = Violin(Ne.Tval(in), i);
                % adjust the violin plot
                misc.violin_setting(v)
                % change colors
                v.ScatterPlot.MarkerFaceColor = col_pos_neg(2,:);
                v.ScatterPlot.MarkerEdgeColor = col_pos_neg(2,:);
                v.ViolinPlot.FaceColor        = col_pos_neg(2,:);
                % save coordinates for plotting lines
                x_tmp(2,:) = v.ScatterPlot.XData;
                y_tmp(2,:) = v.ScatterPlot.YData;
                % plot lines
                arrayfun(@(j) line([x_tmp(1,j), x_tmp(2,j)],[y_tmp(1,j), y_tmp(2,j)], 'Color',...
                    [col(dic(sign([-y_tmp(1,j)+y_tmp(2,j)]/[-x_tmp(1,j)+x_tmp(2,j)])),:), .25],...
                    'LineWidth', .75)...
                    , 1:length(x_tmp))
                clear v x_tmp y_tmp
                i = i + 1;
            end

            ylim([-5,5]);
            set(gca,'Xtick', [1.5,4.5],... 6.5,8.5
                'XTickLabel', ...
                {'OFC', 'vmPFC'},...
                'LineWidth', 2, 'Ytick', [-5:5:5], 'FontSize', 20, 'FontName', 'Arial')
            ylabel('HFB (t-value)')
            xlim([0,6])
            line(xlim(),[0 0], 'LineStyle', '--', 'Color', [.2, .2, .2, .6])

            axis square
            print -dpng -r300 results\Fig1F.png
            print -dsvg results\Fig1F.svg

        end % violin

        function obj = bars(obj, H)
            % bars plot the bar garphs of LMM coefficient
            % Input:
            %      -H   : a string (off/on) determines if the bar graph
            %      layout
            % ------------------------------------------------
            if nargin<2
                H = 'off';
            end
            % mdl bootstrapped model
            Coeff = cat(2,obj.mdl.Coeff)';

            obj = obj.parseCoeff();
            % create the id indices
            id = [
                find(obj.index.OFC & obj.index.Po)  % OFC and positive connotation
                find(obj.index.OFC & obj.index.Ne)  % OFC and negative connotation
                find(~obj.index.OFC & obj.index.Po) % vmPFC and positive connotation
                find(~obj.index.OFC & obj.index.Ne) % vmPFC and negative connotation
                ];
            figure % open a new figure
            c = 0; % set counter
            hold on
            clear b ticks % clear b and ticks
            col = obj.colors([1,2], 2); % get two colors with two reps
            for id_ = 1:length(id) % one: self-coherent minus one: self-incoherent
                c = c+1; % increment the counter by one

                b(c) = bar(c,  quantile(Coeff(:,id(id_,:)),.5), 'Horizontal', H); % plt the actubal bar
                b(c).FaceColor =col(id_,:); % change the bar color to the corresponding colors
                % add the error bar
                stat.errorbar_bootstrap(Coeff, id, id_, c, H)
                if c ==2 % if counter equals 2 additionally increment the counter for adding gaps to bars for visual purposes
                    c = c+1;
                end % end if
            end % end for


            if strcmp(H, 'off')
                xlim([0,6])
                ylim([-.85,1.25])

                set(gca,'Xtick', [1.5,4.5],..., 6.5,8.5
                    'XTickLabel', ...
                    {'OFC', 'vmPFC'},...
                    'LineWidth', 2, 'FontName', 'Arial Nova Cond', 'FontSize', 20)
                set(gca, 'Ytick', sort([0,ylim]))
                axis square
                legend(b(1:2), {'Positive connotation', 'Negative connotation'},...
                    'FontName', 'Arial Nova Cond','Location','northoutside','box', 'off')
            else
                yl =ylim()*1.6;
                % % pbaspect([1,.25,1])
                ylim([0,6])
                xlim([-2,2])
                % % rectangle('Position', [-4.1,  yl(1), .2,  abs(diff(yl))],'FaceColor', 'w', 'EdgeColor', 'none')
                set(gca,'YTick', [1.5,4.5],... 6.5,8.5
                    'YTickLabel',{'OFC','vmPFC'},...
                    'XTick', [-2:2:2],...misc.transormation_scale([-14, -9:3:6],th)
                    'XTickLabel', [-2:2:2],...
                    'LineWidth', 2, 'FontName', 'Arial Nova Cond', 'FontSize', 20)
                % set(gca, 'Ytick', sort([0,ylim]))
                ylabel({'BDI', 'bootstrapped coefficient (%95 CI)'})
                xlabel('HFB (T-value)','FontName', 'Arial Nova Cond', 'FontSize', 18)
                legend(b(1:2), {'Positive connotation', 'Negative connotation'},...
                    'FontName', 'Arial Nova Cond','Location','northoutside','box', 'off')

            end

            if ~exist('results\Fig2F.svg', 'file')
                print -dsvg results\Fig1F2.svg
                print -dpng -r300 results\Fig1F2.png
            else
                print -dsvg results\Fig2A.svg
                print -dpng -r300 results\Fig2A.png
            end
        end % bars

        function obj = parseCoeff(obj)
            % parseCoeff parse the coefficient's name for the given model
            obj.index.OFC = cellfun(@(x) contains(x,'JPAnatomy_OFC'), obj.mdl(1).CoeffName); % anatomical sites
            obj.index.Po  = cellfun(@(x) contains(x,'task_one'), obj.mdl(1).CoeffName); % positive connotation
            obj.index.Ne  = cellfun(@(x) contains(x,'task_minusone'), obj.mdl(1).CoeffName); % negative connotation
        end
        function [POS, NEG, POSXNEG] = subWPoNe(obj)
            % subWPoNe finds individuals who have both positive and
            % negative trials
            % Output:
            %        POS:     a table containing the data for positive
            %                 condition from subjects  with both positive and
            %                 negative
            %        NEG:     a table containing the data for negative
            %                 condition from subjects  with both positive and
            %                 negative
            %        POSXNEG: an intersection of POS and NEG
            % --------------------------------------------

            POS      = []; % initializing the POS table
            NEG      = []; % initializing the NEG table
            for s = unique(obj.preprocT.subj)' % for every unique subjects
                if any(strcmp(obj.preprocT.subj, s{:}) & obj.preprocT.task == 'one') % checks if there is any pos condition for a given subject
                    POS = [POS
                        stat.average_over_rois(stat.average_over_sessions(obj.preprocT(strcmp(obj.preprocT.subj, s{:}) & ...
                        obj.preprocT.task == 'one',:)))];
                end
                if any(strcmp(obj.preprocT.subj, s{:}) & obj.preprocT.task == 'minusone')  % checks if there is any neg condition for a given subject
                    NEG = [NEG
                        stat.average_over_rois(stat.average_over_sessions(obj.preprocT(strcmp(obj.preprocT.subj, s{:}) & ...
                        obj.preprocT.task == 'minusone',:)))];
                end
            end % for every unique subjects
            % average over patients for POS
            pos_s = stat.average_over_subj(POS);
            % remove the enteries with unvalid BDI
            pos_s(isnan(pos_s.BDI),:) = [];
% %             pos_s(pos_s.BDI == 0,:) = [];

            % average over patients for NEG
            neg_s = stat.average_over_subj(NEG);
            % remove the enteries with unvalid BDI
            neg_s(isnan(neg_s.BDI),:) = [];
% %             neg_s(neg_s.BDI == 0,:) = [];
            % find the intersection of two tables
            [~,ip,in] = intersect(pos_s.subj, neg_s.subj);
            % collect the intersection in a cell
            POSXNEG = {pos_s(ip,:), neg_s(in,:)};
        end

        %----- get methods
        function pT = get.preprocT(obj)
            disp('preprocssing the input table...')
            pT = obj.preprocessing(obj.data); % preprcossing the table
            disp('done!')
        end % get method for preprocT
    end %method

    methods(Static)
        function preprocT = preprocessing(Table)
            preprocT = Table;
            % converting all the cell string data to categorical variable
            preprocT.task      = categorical(preprocT.task, {'one','minusone'}, 'Ordinal',true);
            preprocT.JPAnatomy = categorical(preprocT.JPAnatomy);
            % -------- add hemisphere side to data
            preprocT.hemi(preprocT.X>0) = {'r'};
            preprocT.hemi(preprocT.X<0) = {'l'};
            preprocT.hemi = categorical(preprocT.hemi);

            % -------- calculating density across the anatomical sites
            for hemi = ["l", "r"]
                for jp = ["OFC", "MPFC"]
                    preprocT.Density(preprocT .hemi == hemi & preprocT.JPAnatomy == jp ) = ...
                        sum(preprocT.hemi == hemi & preprocT.JPAnatomy == jp );
                end % for sites
            end % for hemis

            preprocT.act  = preprocT.Pval<=.05;
            preprocT.BDIz = (preprocT.BDI - nanmean(preprocT.BDI))./nanstd(preprocT.BDI);
            %  remove the outliers 5
            preprocT(abs(preprocT.Tval)>5, :) =[];
        end % preprocessing

        function col = colors(index, rep)
            % get colors for plotting
            % Input:
            %       -index: a column vector of color indices.
            %       -rep:   a scalar indicates how many times colors should
            %               be repeated
            % Output:
            %       -col:   M X 3 color matrix
            % --------------------------------------------------------------
            arguments % check the input arguments
                index(:,1) double
                rep (1,1) double
            end
            col_resource =  [ 0.1597    0.7772         0
                0.6745         0    0.2329];
            if index> size(col_resource,1) % check if color resource has enough colors
                error('only 6 colors are available, use rep to repeate colors if more are needed!')
            end
            col = col_resource(repmat(index, rep, 1),:);
        end % colors
        function [msr, BDI, b, tval] = OFCvmPFCMEASURE(POSNEG)
            % computes the OFC vmPFC mewasure
            % Input:
            %       -POSNEG:  a 1x2 cell containing the tabular data for
            %                   pos and neg condition HFB as well BDI for same individuals.
            % Output:
            %       - msr:     a vector contaning the OFC vmPFC measure
            %       - BDI:     a vector of BDI scores
            %       - b:       a vector of coeffcients for regressing the
            %                  hemisphere side
            % -----------------------------------------------------
            for i = 1:length(POSNEG)
                [b(:,i), ~, msr_(:,i)] = regress(POSNEG{i}.Tval,...
                    [double(POSNEG{i}.X<0)]); % controlling for hemisphere side
            end
            if all(POSNEG{1}.BDI == POSNEG{2}.BDI) % make sure that the two BDI corresponds together
                BDI = POSNEG{1}.BDI;
            else
                error('mistmach between the BDIs of POS and NEG, something went wrong!')
            end
            tval = msr_;
            msr_ =stat.LMMPoNe.rescale(msr_');% rescalling and sfotmaxing
            msr = msr_(1,:).';
        end % OFCvmPFCMEASURE

        function [curve_score, goodness_score, output_score, pValue] = plot(msr, BDI, tval)
            figure % open a new figure
            hold on
            for ii = 1:length(msr)
                scatter(BDI(ii), msr(ii) , 120, 'kx', 'LineWidth', 2, 'MarkerEdgeAlpha', .75)
            end

            [BDIs, ix] = sort(BDI);
            

            % adjust axis limits
            xlim([0, 40])
            ylim([-1,1])
            % axis ratio
            pbaspect([1,.5,1])
            % draw the plot now
            drawnow()
            % keep the plot and add the regression line
            hold on
            option = fitoptions('poly1');
            % fit a line
            [curve_score, goodness_score, output_score] = fit(msr(ix), BDIs,'poly1', option); 
            [curve_pos, goodness_pos, ~] = fit(tval(ix,1), BDIs,'poly1', option); 
            [curve_neg, goodness_neg, ~] = fit(tval(ix,2), BDIs,'poly1', option); 
            
            % perform anova to get p-value
            % Obtain the coefficient estimates and confidence intervals
            coeff = coeffvalues(curve_score);
            confInt = confint(curve_score);
            % Extract p-value from confidence intervals
            tcritical = tinv(.975, output_score.numobs-output_score.numparam);
            pValue = 2 * (1 - tcdf(abs(coeff./(diff(confInt)/(4*tcritical))),...
                output_score.numobs-output_score.numparam));
            % add the dash line
            col = stat.LMMPoNe.colors([1,2], 2); % get two colors with two reps
            l(1) = dashline(curve_score(-1:.01:1), -1:.01:1, 2, 3, 2, 3, 'LineWidth', 1.25, 'Color', 'k');
                      % add color patches
            yl = ylim();
            p = [0,13.5
                13.5,19.5
                19.5,28.5
                28.5,63];

            color = [0.244, 0.525, 0.310
                0.838, 0.655, 0.509
                0.818, 0.439, 0.380
                0.644, 0.200, 0.298];
            for i= 1:length(p)
                ptc(i) =patch([p(i,:), fliplr(p(i,:))], [yl(1), yl(1), yl(2), yl(2)], color(i,:), 'EdgeColor', 'None', 'FaceAlpha', .25);
            end
           

            set(gca, 'FontName', 'Arial Nova Cond', 'FontSize', 18, 'LineWidth', 2);
            ylabel('OFC&vmPFC Score')
            xlabel('BDI')

            yyaxis right 
            ax= gca();
            ax.YAxis(2).Color = 'k';
            ylabel('HFB (t-value)')
            hold on
           if goodness_pos.adjrsquare >0
             l(2) = dashline(curve_pos(-5:.01:5), -5:.01:5, 1, 4, 1, 4, 'LineWidth', 1.5, 'Color', col(1,:));
           else
              l(2) = dashline(linspace(0, 40, length((-5:.01:5))), 0.*(-5:.01:5), 1, 4, 1, 4, 'LineWidth', 1.5, 'Color', col(1,:));
           end
           if goodness_neg.adjrsquare >0
            l(3) = dashline(curve_neg(-5:.01:5), -5:.01:5, 1, 4, 1, 4, 'LineWidth', 1.5, 'Color', col(2,:));
           else
            l(3) = dashline(linspace(0, 40, length((-5:.01:5))), 0.*(-5:.01:5), 1, 4, 1, 4, 'LineWidth', 1.5, 'Color', col(2,:));
           end
            ylim([-5,5])

            print -dpng -r300 results\Fig2B.png
            print -dsvg results\Fig2B.svg
        end % plot
        function y = softmax(x)
            y = exp(x) ./ sum(exp(x));
        end % softmax
        function y = rescale(x)
            y = 2.*(stat.LMMPoNe.softmax(x)-.5);

        end % rescale 
    end % methods(Static)
end % class
% $END

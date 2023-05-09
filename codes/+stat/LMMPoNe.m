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

 

        function obj = bars(obj)
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
                b(c) = bar(c,  quantile(Coeff(:,id(id_,:)),.5)); % plt the actubal bar
                b(c).FaceColor =col(id_,:); % change the bar color to the corresponding colors
                % add the error bar
                errorbar(c, quantile(Coeff(:,id(id_,:)),.5), ... median
                    -quantile(Coeff(:,id(id_,:)),.5)+quantile(Coeff(:,id(id_,:)),.025), ... low boundry of CI
                    quantile(Coeff(:,id(id_,:)),.975)-quantile(Coeff(:,id(id_,:)),.5), ... high boundry of CI
                    'Color', 'k', 'Linewidth', 1, 'CapSize', 0)

                if c ==2 % if counter equals 2 additionally increment the counter for adding gaps to bars for visual purposes
                    c = c+1;
                end % end if
            end % end for
            pbaspect([.75,1,1])
            xlim([0,6])

            set(gca,'Xtick', [1.5,4.5],...
                'XTickLabel', ...
                {'OFC','vmPFC',},...
                'LineWidth', 2, 'FontName', 'Arial Nova Cond', 'FontSize', 20)
            set(gca, 'Ytick', sort([0,ylim]))
            legend(b(1:2), {'Positive connotation', 'Negative connotation'},...
                'FontName', 'Arial Nova Cond','Location','eastoutside','box', 'off')

            print -dsvg results\Fig1F.svg
            print -dpng -r300 results\Fig1F.png
        end % bars

        function obj = parseCoeff(obj)
            % parseCoeff parse the coefficient's name for the given model
            obj.index.OFC = cellfun(@(x) contains(x,'JPAnatomy_OFC'), obj.mdl(1).CoeffName); % anatomical sites
            obj.index.Po  = cellfun(@(x) contains(x,'task_one'), obj.mdl(1).CoeffName); % positive connotation
            obj.index.Ne  = cellfun(@(x) contains(x,'task_minusone'), obj.mdl(1).CoeffName); % negative connotation
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
    end % methods(Static)
end % class


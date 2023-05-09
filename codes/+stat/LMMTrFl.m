classdef LMMTrFl < stat.LMM
    %   CLASS NAME: LMMSTrFl
    %
    %   Purpose: Linear mixed-effect model for comparing self-coherence in
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
    %   Date: 05/08/2023
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        index  % a structure containing the indices for each coefficient
    end
    properties(Dependent)
        preprocT % preprocessed data for True(self-coherent) and False (self-incoherent) LMM
    end
    methods
        function obj = bars(obj)
             % mdl bootstrapped model
            Coeff = cat(2,obj.mdl.Coeff)';
            
             obj = obj.parseCoeff();
            % create the id indices
            id = [
                find(obj.index.OFC & obj.index.SCoh)  % OFC and self-coherent
                find(obj.index.OFC & obj.index.SIcoh)  % OFC and self-incoherent
                find(~obj.index.OFC & obj.index.SCoh) % vmPFC and self-coherent
                find(~obj.index.OFC & obj.index.SIcoh) % vmPFC and self-incoherent
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
                stat.errorbar_bootstrap(Coeff, id, id_, c, 'off'); % adding CI to the bar plots  
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
            legend(b(1:2), {'self-coherent', 'self-incoherent'},...
                'FontName', 'Arial Nova Cond','Location','eastoutside','box', 'off')

            print -dsvg results\FigS2.svg
            print -dpng -r300 results\FigS2.png
        end % bars

        function obj = parseCoeff(obj)
            % parseCoeff parse the coefficient's name for the given model
            obj.index.OFC = cellfun(@(x) contains(x,'JPAnatomy_OFC'), obj.mdl(1).CoeffName); % anatomical sites
            obj.index.SCoh = cellfun(@(x) contains(x,'task_true'), obj.mdl(1).CoeffName); % self-coherent
            obj.index.SIcoh = cellfun(@(x) contains(x,'task_false'), obj.mdl(1).CoeffName); % self-incoherent
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
            preprocT.task      = categorical(preprocT.task, {'true','false'}, 'Ordinal',true);
            preprocT.task  = reordercats(preprocT.task, {'true', 'false'});
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
            col_resource =  [.5, .5, .65
                                .65, .5, .5];
            if index> size(col_resource,1) % check if color resource has enough colors
                error('only 6 colors are available, use rep to repeate colors if more are needed!')
            end
            col = col_resource(repmat(index, rep, 1),:);
        end % colors
    end % methods(Static)
end % class


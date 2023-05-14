classdef supplmentary
    %   CLASS NAME: supplmentary
    %
    %   Purpose: provides for supplemntary statistics
    %   This script is part of the codes for reproducing stats reported in
    %   "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
    %
    %   Properties:
    %               - LMM:      a class of LMM
    %   Methods:
    %
    %   Copyright (C)  Behzad Iravani, department of neurology and neurological
    %   sciences, Stanford University. May 2023.
    %
    %   Author: Behzad Iravani
    %   behzadiravani@gmail.com
    %   Contact: behzadiravani@gmail.com
    %   Date: 05/10/2023
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        LMM
    end % properties
    methods
        function obj = supplmentary(LMM)
            obj.LMM = LMM;
        end
        function histg(obj, x, g)
            levels = unique(obj.LMM.data.(g));
            nplot = length(levels);
            for iplot = 1:nplot
                subplot(nplot,1,iplot)
                histogram(obj.LMM.data.(x)(strcmp(obj.LMM.data.(g), levels{iplot})), 5, FaceColor=obj.LMM.colors(iplot,1)) % positive RT

                axis square
                box off
                set(gca, 'FontSize', 16, ...
                    'FontName', 'Arial Nova Cond', 'Xtick', 1:5)
            end % for iplot
            ylabel('Counts')
            xlabel(x)
        end % histg

        function scatterg(obj, x, y, g)

            figure

            hold on
            levels = unique(obj.LMM.data.(g));
            nplot = length(levels);
            for iplot = 1:nplot

                s(iplot) = scatter(...
                    obj.LMM.data.(x)(strcmp(obj.LMM.data.(g), levels{iplot})),...
                    obj.LMM.data.(y)(strcmp(obj.LMM.data.(g), levels{iplot})),...
                    100, MarkerFaceColor=obj.LMM.colors(iplot,1),...
                    MarkerEdgeColor=obj.LMM.colors(iplot,1));

                set(gca, 'FontSize', 16, ...
                    'FontName', 'Arial Nova Cond', 'LineWidth', 2)
            end % for iplot
            xlabel(x)
            ylabel(y)
        end % scatterg
        function mdl = stat(obj, model)
            % stat computes supplmentary statistics using data in LMM
            % object.
            % --------------------------------
            arguments
                obj (1,1)
                model string
            end
            mdl = fitlme(obj.LMM.preprocT, model);


        end % stats

        function scat_his(obj)

            sch = scatterhist(obj.LMM.preprocT.BDI(obj.LMM.preprocT.BDI~=0),...
                obj.LMM.preprocT.Tval(obj.LMM.preprocT.BDI~=0),'Group',...
                obj.LMM.preprocT.hemi(obj.LMM.preprocT.BDI~=0),'Marker','sd', 'Kernel','on', 'Bandwidth', [3,3;.3,.3]);

            axes(sch(1))
            box off
            set(gca, 'LineWidth',2, 'FontName', 'Arial Nova Cond', 'FontSize', 18, 'TickLength', [.02,.02])
            xlabel('BDI')
            ylabel('HFB (t-value)')

            legend({'r-hemi','l-hemi'}, 'Box', 'off', 'FontName', 'Arial Nova Cond')
            axis square
            xlim([0, 40])
            ylim([-4,4])
            print -dpng -r300 results\supp--BDI--hemi.png
            print -dsvg results\supp--BDI--hemi.svg
        end % scat_hist
        function report(obj)

            Coeff = cat(2,obj.LMM.mdl.Coeff)';

            id = logical(dummyvar(1:size(Coeff,2)));
            for icolumn=1:size(Coeff,2)
                m(icolumn) = mean(Coeff(:,icolumn));
                s(icolumn)= std(Coeff(:,icolumn));

                [lowerbound(icolumn), upperbound(icolumn)] =stat.errorbar_bootstrap(Coeff, id, icolumn, icolumn, 'off');
                CI(:,(icolumn)) =   [m(icolumn) - lowerbound(icolumn)
                    upperbound(icolumn) - m(icolumn)];
                close
                p(icolumn) =   (1+sum(sign(mean(Coeff(:,id(icolumn,:))))*Coeff(:,id(icolumn,:))<=0))./(length(Coeff)+1);
                
                fprintf('%s: %1.2f +/- %1.2f ci = [%1.2f, %1.2f], p = %.3g\n', obj.LMM.mdl(1).CoeffName{id(icolumn,:)}, ...
                    m(icolumn), ...
                    m(icolumn), ...
                    CI(1,icolumn), CI(2,icolumn),...
                    p(icolumn)...
                );
            end


        end % bars
    end % methods
end %
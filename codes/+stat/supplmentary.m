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
        function mdl = stat(obj, model, dummy, avgoverroi, avgovertask, avgoversubj)
            % stat computes supplmentary statistics using data in LMM
            % object.
            % --------------------------------

            if nargin<3
                dummy = 'reference';
                avgoverroi= false;
                avgovertask = false;
                avgoversubj = false;
            end

            tmpT = obj.LMM.preprocT;
            if avgoverroi
                tmpT = stat.average_over_rois(tmpT);
            end
            if avgovertask
                icond = 0;
                for conds = unique(tmpT.task)'
                    icond = icond +1;
                    tmpT_{icond} = stat.average_over_task(tmpT(tmpT.task == conds,:));
                    tmpT_{icond}.task = repmat(conds, height(tmpT_{icond}),1);
                end
            end
            if avgoversubj
                for i = 1:length(tmpT_)
                    tmpT_{i} = stat.average_over_subj(tmpT_{i});
                end
                tmpT = cat(1, tmpT_{:});
            end

            mdl = fitlme(tmpT, model, 'DummyVarCoding', dummy);


        end % stats

        function scat_his(obj, group)
            warning off
            switch group
                case "hemi"
                    sch = scatterhist(obj.LMM.preprocT.BDI,...
                        obj.LMM.preprocT.Tval,'Group',...
                        obj.LMM.preprocT.hemi,'Marker','sd', 'Kernel','on', 'Bandwidth', [3,3;.3,.3]);

                    %             [H.t, pValue.t, KSstatistic.t] = kstest2(obj.LMM.preprocT.Tval( obj.LMM.preprocT.hemi == 'l'),...
                    %                 obj.LMM.preprocT.Tval( obj.LMM.preprocT.hemi == 'r'), 'Tail','unequal');
                    %
                    %
                    %             [H.bdi, pValue.bdi, KSstatistic.bdi] = kstest2(obj.LMM.preprocT.BDI( obj.LMM.preprocT.hemi == 'l'),...
                    %                 obj.LMM.preprocT.BDI( obj.LMM.preprocT.hemi == 'r'), 'Tail','unequal');

                    [p.BDI,anovatab.BDI,stats.BDI] = anova1(obj.LMM.preprocT.BDI, obj.LMM.preprocT.hemi);
                    [p.T,anovatab.T,stats.T] = anova1(obj.LMM.preprocT.Tval, obj.LMM.preprocT.hemi);

                    axes(sch(1))
                    box off
                    set(gca, 'LineWidth',2, 'FontName', 'Arial Nova Cond', 'FontSize', 18, 'TickLength', [.02,.02])
                    xlabel('BDI')
                    ylabel('HFB (t-value)')

                    legend({'r-hemi','l-hemi'}, 'Box', 'off', 'FontName', 'Arial Nova Cond')
                    axis square
                    xlim([0, 40])
                    ylim([-4,4])
                case "valence"
                    sch = scatterhist(obj.LMM.preprocT.BDI,...
                        obj.LMM.preprocT.Tval,'Group',...
                        obj.LMM.preprocT.task,'Marker','sd', 'Color', stat.LMMPoNe.colors([1,2],1),...
                        'Kernel','on', 'Bandwidth', [3,3;.3,.3]);
                    

                    vdat = stat.average_over_sessions(obj.LMM.preprocT);
                    Po =  stat.average_over_task(vdat(strcmp(vdat.task,'one'),:)); % Positive data
                    Ne =  stat.average_over_task(vdat(strcmp(vdat.task,'minusone'),:)); % Negative data;
                    % ......
                    [~, ip, in] = intersect(cellfun(@(x,y) sprintf('%s:%s',x,y), Po.subj, Po.chan, 'UniformOutput', false),...
                        cellfun(@(x,y) sprintf('%s:%s',x,y), Ne.subj, Ne.chan, 'UniformOutput', false)); % include those that have value in both tables

                    Po = Po(ip, :);
                    Ne = Ne(in,:);
                    dataTable = table();
                    dataTable = [Po; Ne];
                    dataTable.Cond = [repmat({'pos'}, height(Po),1); repmat({'neg'}, height(Po),1)];
%                     dataTable.T1 = Po.Tval;
%                     dataTable.T2 = Ne.Tval;
%                     dataTable.BDI = Po.BDI;
%                     dataTable.subject = Po.subj;
                    
                    % Perform repeated measures ANOVA
                    
                    lmModel = fitlme(dataTable, 'Tval~ 1 + BDI*Cond + (1|subj)');
                    [T] = anova(lmModel);

                    % % %                     [p.BDI,anovatab.BDI,stats.BDI] = anovan(obj.LMM.preprocT, obj.LMM.preprocT.task);
                  
                   
                    drawnow()
                    axes(sch(1))
                    box off
                    set(gca, 'LineWidth',2, 'FontName', 'Arial Nova Cond', 'FontSize', 18, 'TickLength', [.02,.02])
                    xlabel('BDI')
                    ylabel('HFB (t-value)')
                   
                    legend({'Postive Con.','Negative Con'}, 'Box', 'off', 'FontName', 'Arial Nova Cond')
                    axis square
                    xlim([0, 40])
                    ylim([-4,4])

                    bdi = 1:35;
                    ax = axes('Position', sch(1).Position, 'InnerPosition', sch(1).InnerPosition);
                    axis square off
                    uistack(ax,"top")
                   
                    line(bdi, lmModel.Coefficients.Estimate(1) + lmModel.Coefficients.Estimate(2)*bdi,...
                        'Color', stat.LMMPoNe.colors([1],1), 'LineWidth', 2) % positive
                    line(bdi, lmModel.Coefficients.Estimate(1) + lmModel.Coefficients.Estimate(3)+...
                        (lmModel.Coefficients.Estimate(2) + lmModel.Coefficients.Estimate(4))*bdi,...
                        'Color', stat.LMMPoNe.colors([2],1), 'LineWidth', 2) % positive
                        delete(sch(2));
            jsontext = jsonencode(struct(T), 'PrettyPrint',true);
            fid = fopen('results\figS6_stat.json', 'w');
            fprintf(fid, jsontext);
            fclose(fid);
            end 
            
           
%         

            drawnow()
            print(strcat("results\supp--BDI--",group,".png"), '-dpng', '-r300')
            print(strcat("results\supp--BDI--",group,".svg"), '-dsvg')
            warning on

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
%-*- UTF8 -*-
classdef stat_report
    % CLASS NAME: stat_report
    %
    % Purpose: stat_report provides methods for reporting the demographic,
    % electrodes' coverage and etc. for the manuscript
    % "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
    %
    % Properties:
    %   - data: a table containing the summary of data
    %
    % Methods:
    %   - report(what): report the info based on the input (what),
    %   and behvaioral .json file see below for more info.
    %
    % Author: Behzad Iravani
    % Contact: behzadiravani@gmail.com
    % Date: 05/01/2023
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        data % a table containig data -- see bellow for more information
        jbpath  %  path to *.json behavioral data
        jepath  %  path to *.json electrode type file
    end

    properties(Dependent)
        BHV         % structure containing the behavioral data
        ECoGSEEG    % structure containing electrode types
        indexSelf    
        indexNotSelf
    end

    methods
        function obj = stat_report(data, jbpath, jepath)
            % The constructor method creates the instanse of the class.
            % Input:
            %       data: a table containig data
            %
            %               subj      chan      task        X         Y          Z       dof    responseTr    Loc_Tval    Loc_Pval         avg             time           RT      varRT     JPAnatomy      Tval        Pval
            %     _____    _______    _____    _______    ______    _______    ___    __________    ________    ________    _____________    _____________    ______    ______    _________    _________    _______
            %
            %     "S20"    "AFS10"    "EP" (Self-episodic)     -16.596    68.162    -1.6345    24         11         -1.0283    0.073385    1×2751 double    1×2751 double    4.2617    4.3096     "MPFC"        -1.1689    0.12697
            %     "S20"    "AFS10"    "MTH"    -16.596    68.162    -1.6345    36         18             NaN         NaN    1×2751 double    1×2751 double    5.7958    15.186     "MPFC"       -0.96219    0.17596
            %     "S20"    "AFS10"    "OTH"    -16.596    68.162    -1.6345    24         20          -1.234    0.014197    1×2751 double    1×2751 double     5.418    3.2609     "MPFC"       -0.59045    0.28754
            %
            %       jpath:  a string contais the path to behvioral json file

            % **************************************************************************************************************************************************************************************************************
            obj.data   = data;
            obj.jbpath = jbpath;
            obj.jepath = jepath;
        end % stat_report: constructor method

        function  out = report(obj, what)
            % REPORT generates a short information for given argument (what).
            % Input:
            %       what: a string determining the type of report to be
            %              genrated:
            %              -"num_indiv": reports the neumber of unique
            %               individual in the dataset.
            %
            %              -"number_total_elec": report the number of total
            %               unique electrodes in the dataset.
            %
            %              -"number_trials": computes summary statistics of
            %               number of trials
            %
            %              -"number_true_false"
            %
            %              -"reaction_time"
            %
            %              -"veridicality"
            %
            %              -"ECoGSEEG"
            % ------------------------------------------------------------------------------
            switch what
                case "num_indiv"
                    subj = unique(obj.data.subj); % Get the unique subject IDs from the data
                    fprintf('The total number of pt is: %d\n\r ',...
                        length(subj)) % Print the total number of patients in the data to the console

                case "number_total_elec"
                    elec = numbertotalelec(obj); % Get the total number of electrodes for each patient
                    fprintf('The total number of elec is: %d, in total patients %d\n mean (std) # elec: %1.2f(%1.2f), range =  [%1.0f,%1.0f]\n\r ',...
                        sum(elec), numel(elec), mean(elec), std(elec),...
                        min(elec), max(elec)); % Compute and print summary statistics for the number of electrodes

                case "number_trials" % computes summary statistics of number of trials
                    out = cellfun(@(x)[mean(x), std(x)],{obj.BHV.number_ep,...
                        obj.BHV.number_sj ...
                        obj.BHV.number_mth},...
                        UniformOutput = false );% Compute the mean and standard deviation of the number of trials for each behavioral measure

                    cellfun(@(x,y) fprintf('%s # trails: mean (std): %1.0f (%1.1f)\n', x, y), {'EP', 'SJ','MTH'}, out);% Print the mean and standard deviation of the number of trials for each behavioral measure to the console
                case "number_true_false"
                    out.true = cellfun(@(x) round([mean(x), std(x), min(x), max(x)]),{obj.BHV.number_ep_true,...
                        obj.BHV.number_sj_true ...
                        obj.BHV.number_mth_true},...
                        UniformOutput = false );% Compute the mean and standard deviation of the number of trials replied with true

                    out.false = cellfun(@(x) round([mean(x), std(x), min(x), max(x)]),{obj.BHV.number_ep_false,...
                        obj.BHV.number_sj_false ...
                        obj.BHV.number_mth_false},...
                        UniformOutput = false );% Compute the mean and standard deviation of the number of trials replied with false

                    cellfun(@(x,y) fprintf('%s # trails replied with true: mean (std): %1.0f (%1.0f), range = [%1.0f,%1.0f]\n', x, y), {'EP true', 'EP false'}, [out.true(1), out.false(1)]);% Print the mean and standard deviation of the number of trials replied with true and false for each (EP as refered to SE (self-episodic) in the paper, not be confused with SE (Semantic) in the BHV data) measure to the console
                    cellfun(@(x,y) fprintf('%s # trails replied with true: mean (std): %1.0f (%1.0f), range = [%1.0f,%1.0f]\n', x, y), {'SJ true', 'SJ false'}, [out.true(2), out.false(2)]);% Print the mean and standard deviation of the number of trials replied with true and false for each SJ measure to the console
                    cellfun(@(x,y) fprintf('%s # trails replied with true: mean (std): %1.0f (%1.0f), range = [%1.0f,%1.0f]\n', x, y), {'MTH true', 'MTH false'}, [out.true(3), out.false(3)]);% Print the mean and standard deviation of the number of trials replied with true and false for each MTH measure to the console
                case "reaction_time"
                    out.true = cellfun(@(x) round([nanmean(x), nanstd(x), min(x), max(x)],2),{obj.BHV.number_ep_rt_t,...
                        obj.BHV.number_sj_rt_t ...
                        obj.BHV.number_mth_rt_t},...
                        UniformOutput = false );% Compute the mean and standard deviation of the RT replied with true

                    out.false = cellfun(@(x) round([nanmean(x), nanstd(x), min(x), max(x)],2),{obj.BHV.number_ep_rt_f,...
                        obj.BHV.number_sj_rt_f ...
                        obj.BHV.number_mth_rt_f},...
                        UniformOutput = false );% Compute the mean and standard deviation of the RT of trials replied with false

                    cellfun(@(x,y) fprintf('%s RT replied with true: mean (std): %1.2f (%1.2f), range = [%1.2f,%1.2f]\n', x, y), {'EP true', 'EP false'}, [out.true(1), out.false(1)]);% Print the mean and standard deviation of RT replied with true and false for each SE(EP) measure to the console
                    cellfun(@(x,y) fprintf('%s RT replied with true: mean (std): %1.2f (%1.2f), range = [%1.2f,%1.2f]\n', x, y), {'SJ true', 'SJ false'}, [out.true(2), out.false(2)]);% Print the mean and standard deviation of RT replied with true and false for each SJ measure to the console
                    cellfun(@(x,y) fprintf('%s RT replied with true: mean (std): %1.2f (%1.2f), range = [%1.2f,%1.2f]\n', x, y), {'MTH true', 'MTH false'}, [out.true(3), out.false(3)]);% Print the mean and standard deviation of RT replied with true and false for each MTH measure to the console

                    % prepare data for repeated measure ANOVA
                    warning off
                    datt = table([], {}, {},'VariableNames',{'rt1','task', 'TF'});
                    for cond1 = ["number_ep_rt_t","number_sj_rt_t","number_mth_rt_t"]
                        H = height(datt);
                        datt.rt1(H+1:H+length(obj.BHV.(cond1)))   = obj.BHV.(cond1);
                        datt.task(H+1:H+length(obj.BHV.(cond1))) = ...
                            repmat(....
                            {char(regexp(cond1, '(?<=_)\w+(?=_rt)','match'))},...
                            length(obj.BHV.(cond1)),1);
                        datt.TF(H+1:H+length(obj.BHV.(cond1))) = ...
                            repmat(....
                            {char(regexp(cond1, "(?<=_rt_)\w+",'match'))},...
                            length(obj.BHV.(cond1)),1);
                    end % cond 1

                    datf = table([], {}, {},'VariableNames',{'rt2','task', 'TF'});
                    for cond1 = ["number_ep_rt_f", "number_sj_rt_f", "number_mth_rt_f"]
                        H = height(datf);
                        datf.rt2(H+1:H+length(obj.BHV.(cond1)))   = obj.BHV.(cond1);
                        datf.task(H+1:H+length(obj.BHV.(cond1))) = ...
                            repmat(....
                            {char(regexp(cond1, '(?<=_)\w+(?=_rt)','match'))},...
                            length(obj.BHV.(cond1)),1);
                        datf.TF(H+1:H+length(obj.BHV.(cond1))) = ...
                            repmat(....
                            {char(regexp(cond1, "(?<=_rt_)\w+",'match'))},...
                            length(obj.BHV.(cond1)),1);
                    end % cond 1
                    warning on
                    task = char(datt.task');

                    t = table(task, datt.rt1, datf.rt2, ...
                        'VariableNames',{'task','rt1','rt2'});
                    rm = fitrm(t, ' rt1-rt2 ~ task', 'WithinDesign', 1:2, 'WithinModel','orthogonalcontrasts');
                    [A.tbl,A.A,A.C,A.D] = ranova(rm);
                    % compute rmANOVA statistics for RT and self
                    fid = fopen('results\reaction_time.json', 'w');
                    jsontext = jsonencode(A, 'PrettyPrint', true);
                    fprintf(fid, jsontext);
                    fclose(fid);

                    % Self-referntial model
                    rm = fitrm(t(~strcmp(cellstr(t.task), 'mth'),:), ' rt1-rt2 ~ task', 'WithinDesign', 1:2, 'WithinModel','orthogonalcontrasts');
                    [A.tbl,A.A,A.C,A.D] = ranova(rm);
                    % compute rmANOVA statistics for RT and self
                    fid = fopen('results\reaction_time_SelfStats.json', 'w');
                    jsontext = jsonencode(A, 'PrettyPrint', true);
                    fprintf(fid, jsontext);
                    fclose(fid);
                    % math
                    rm = fitrm(t(strcmp(cellstr(t.task), 'mth'),:), ' rt1-rt2 ~ 1', 'WithinDesign', 1:2, 'WithinModel','orthogonalcontrasts');
                    [A.tbl,A.A,A.C,A.D] = ranova(rm);
                    % compute rmANOVA statistics for RT and self
                    fid = fopen('results\reaction_time_MthStats.json', 'w');
                    jsontext = jsonencode(A, 'PrettyPrint', true);
                    fprintf(fid, jsontext);
                    fclose(fid);

                case "veridicality" % computes summary statistics of number of veridicality

                    out.true = cellfun(@(x) round([nanmean(x), nanstd(x), min(x), max(x)],2),{obj.BHV.verdicallity_SE_true,... EP, no veridcality could be determined for SJ
                        obj.BHV.verdicallity_MTH_true},...
                        UniformOutput = false );% Compute the mean and standard deviation of response veridicality responded with true for each behavioral measure

                    out.false = cellfun(@(x) round([nanmean(x), nanstd(x), min(x), max(x)],2),{obj.BHV.verdicallity_SE_false,... EP, no veridcality could be determined for SJ
                        obj.BHV.verdicallity_MTH_false},...
                        UniformOutput = false );% Compute the mean and standard deviation of response veridicality responded with false for each behavioral measure

                    cellfun(@(x,y) fprintf('%s veridicality replied with true: mean (std): %1.2f (%1.2f), range = [%1.2f,%1.2f]\n', x, y), {'EP true', 'EP false'}, [out.true(1), out.false(1)]);% Print the mean and standard deviation of veridicality replied with true and false for each SE(EP) measure to the console
                    cellfun(@(x,y) fprintf('%s veridicality replied with true: mean (std): %1.2f (%1.2f), range = [%1.2f,%1.2f]\n', x, y), {'MTH true', 'MTH false'}, [out.true(2), out.false(2)]);% Print the mean and standard deviation of veridicality replied with true and false for each MTH measure to the console
                case "ECoGSEEG"
                    % Compute the mean and standard deviation of the number of ECoG and SEEG electrodes, as well as the number of OFC and vmPFC electrodes
                    cellfun(@(s,e) fprintf('%s -- electype: %s\n', s, e), {obj.ECoGSEEG.subj}, {obj.ECoGSEEG.elec_type}); % write the subject id and elec type to the console

                    fprintf('ECOG = %1.0f +/- %1.0f, [%1.0f, %1.0f]\n', ...
                        mean([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'ECOG').NECOG]), ...
                        std([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'ECOG').NECOG ]), ...
                        min([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'ECOG').NECOG ]), ...
                        max([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'ECOG').NECOG ]))
                    fprintf('OFC = %1.2f +/- %1.2f\n', ...
                        mean([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'ECOG').nOFC ]), ...
                        std([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'ECOG').nOFC ]))
                    fprintf('MPFC = %1.2f +/- %1.2f\n', ...
                        mean([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'ECOG').nMPFC]), ...
                        std([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'ECOG').nMPFC]))
                    fprintf('SEEG = %1.0f +/- %1.0f, [%1.0f, %1.0f]\n', ...
                        mean([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'SEEG').NSEEG]), ...
                        std([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'SEEG').NSEEG]), ...
                        min([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'SEEG').NSEEG]), ...
                        max([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'SEEG').NSEEG]))
                    fprintf('OFC = %1.2f +/- %1.2f\n', ...
                        mean([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'SEEG').nOFC]), ...
                        std([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'SEEG').nOFC]))
                    fprintf('MPFC = %1.2f +/- %1.2f\n', ...
                        mean([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'SEEG').nMPFC]), ...
                        std([obj.ECoGSEEG(categorical({obj.ECoGSEEG.elec_type}) == 'SEEG').nMPFC]))
                case "active_total" % computes number of self- and math
                    c = 0;
                    out = {};
                    T = stat.average_over_task(obj.data(strcmp(obj.data.task, 'EP') | strcmp(obj.data.task, 'SJ'),:));
                    index_self = obj.indexSelf;
                    index_not  = obj.indexNotSelf;
                    for s = unique(T.subj)'
                        c = c+1;
                        ch_c = 0;
                        out{c}.self = mean(index_self(categorical(T.subj) == s{:}));
                        index_math = (T.Loc_Pval<.05 & T.Loc_Tval <0)&~index_self;
                        out{c}.math = mean(index_math(categorical(T.subj) == s{:}));
                        out{c}.not = mean(index_not(categorical(T.subj) == s{:}));
                    end
                case "number_regional"

                    out = struct();
                    i_s = 0;
                    for s = unique(obj.data.subj)'
                        i_s = i_s + 1;
                        for jp = unique(obj.data.JPAnatomy)'

                            out.(['l', jp{:}])(i_s) = length(unique(obj.data.chan(strcmp(obj.data.subj, s{:})...
                                & strcmp(obj.data.JPAnatomy, jp{:}) &....
                                obj.data.X<0)));
                            out.(['r', jp{:}])(i_s) = length(unique(obj.data.chan(strcmp(obj.data.subj, s{:})...
                                & strcmp(obj.data.JPAnatomy, jp{:}) &....
                                obj.data.X>0)));
                        end
                    end
                case "active_regions"
                    out = [];
                    for hemi =["left", "right"]
                        for views = ["medial", "ventral"]
                            out.(hemi).(views) = obj.regionals_numbers(hemi,views);
                        end % views
                    end % hemi

                case "percentage"
                    c = 0;
                    out = struct();
                    
                    T = stat.average_over_task(obj.data(strcmp(obj.data.task, 'EP') | strcmp(obj.data.task, 'SJ'),:));
                    index_self = obj.indexSelf;
                    for s = unique(obj.data.subj)'
                        c = c+1;
                        out.nEP{c}  = pre_task_selfact('EP',index_self);
                        out.nSJ{c}  = pre_task_selfact('SJ',index_self);
                        out.both{c} = pre_both_selfact(T, index_self);
                        %                         end
                    end % for

                    case "percentageAll"
                    c = 0;
                   
                    for s = unique(obj.data.subj)'
                        c = c+1;
                        out{c} = pre_all_selfact(obj.indexSelf);
                        %                         end
                    end % for

                     case "percentageAll&"
                    c = 0;
                   
                    for s = unique(obj.data.subj)'
                        c = c+1;
                        out{c} = pre_allAnd_selfact(obj.indexSelf);
                        %                         end
                    end % for
                    
                case "sumTrueFalse"
                    out.true = obj.BHV.sum_number_ep_true  + obj.BHV.sum_number_sj_true;
                    out.false = obj.BHV.sum_number_ep_false  + obj.BHV.sum_number_sj_false;
            end  %switch
            if  exist('subj', 'var')
                out = subj;

            end
            % ------------------------------------------------
            function prec = pre_task_selfact(task, index_self)
               
                index_task = categorical(obj.data.subj) == s{:} &... given subject
                    any(categorical(obj.data.chan) ==  categorical(index_self)',2) &... self-channel
                    strcmp(obj.data.task, task) & ... given task
                    obj.data.Tval>0 &... above baseline
                    obj.data.Pval<.05;

                index_self_subj = categorical(obj.data.subj) == s{:} &... given subject
                    ...strcmp(obj.data.task, task) & ... given task
                    any(categorical(obj.data.chan) ==  categorical(index_self)',2);
                if nansum(index_self_subj)>0
                    prec = nansum(index_task)./nansum(index_self_subj);
                else
                    prec = nan;
                end % if
            end % pre_task_selfact

            function prec = pre_both_selfact(T, index_self)
                index_task = categorical(T.subj) == s{:} &... given subject
                    index_self & ... self-channel
                    T.Tval>0 &... above baseline
                    T.Pval<.05;


                index_self_subj = categorical(T.subj) == s{:} &... given subject
                    index_self;

                if nansum(index_self_subj)>0
                    prec = nansum(index_task)./nansum(index_self_subj);
                else
                    prec = nan;
                end % if
            end % pre_task_selfact
                

            function prec = pre_all_selfact(index_self)
                SE =  stat.average_over_task(obj.data(strcmp(obj.data.task, 'EP'),:));
                
                index_SE =  categorical(SE.subj) == s{:} &... given subject
                    index_self &...
                    SE.Tval>0 &... above baseline
                    SE.Pval<.05;
                
                SJ =  stat.average_over_task(obj.data(strcmp(obj.data.task, 'SJ'),:));
                
                index_SJ =  categorical(SJ.subj) == s{:} &... given subject
                    index_self &...    
                    SJ.Tval>0 &... above baseline
                    SJ.Pval<.05;
                
                Mth =  stat.average_over_task(obj.data(strcmp(obj.data.task, 'MTH'),:));
                
                index_Math =  categorical(Mth.subj) == s{:} &... given subject
                    index_self &...    
                     Mth.Tval>0 &... above baseline
                    Mth.Pval<.05;

%                 index_self_subj = categorical(T.subj) == s{:} &... given subject
%                     index_self;
                if ~((sum(strcmp(SE.subj, s{:})) == sum(strcmp(SJ.subj, s{:}))) &&...
                        (sum(strcmp(SE.subj, s{:})) == sum(strcmp(Mth.subj, s{:}))) &&...
                        (sum(strcmp(SJ.subj, s{:})) == sum(strcmp(Mth.subj, s{:}))))
                        error('Subjects number are not matching across tasks!')
                end
                
                TAll= (index_self & index_SE) | (index_self & index_SJ) |(index_Math &  index_self);
                index_subj  = strcmp(SE.subj, s{:}) &  index_self;

                if nansum(index_subj)>0
                    prec = nansum(TAll)./nansum(index_subj);
                else
                    prec = nan;
                end % if
            end % pre_all_selfact

            function prec = pre_allAnd_selfact(index_self)
                SE =  stat.average_over_task(obj.data(strcmp(obj.data.task, 'EP'),:));
                
                index_SE =  categorical(SE.subj) == s{:} &... given subject
                    SE.Tval>0 &... above baseline
                    SE.Pval<.05;
                
                SJ =  stat.average_over_task(obj.data(strcmp(obj.data.task, 'SJ'),:));
                
                index_SJ =  categorical(SJ.subj) == s{:} &... given subject   
                    SJ.Tval>0 &... above baseline
                    SJ.Pval<.05;
                
                Mth =  stat.average_over_task(obj.data(strcmp(obj.data.task, 'MTH'),:));
                
                index_Math =  categorical(Mth.subj) == s{:} &... given subject  
                    Mth.Tval>0 &... above baseline
                    Mth.Pval<.05;

%                 index_self_subj = categorical(T.subj) == s{:} &... given subject
%                     index_self;
                if ~((sum(strcmp(SE.subj, s{:})) == sum(strcmp(SJ.subj, s{:}))) &&...
                        (sum(strcmp(SE.subj, s{:})) == sum(strcmp(Mth.subj, s{:}))) &&...
                        (sum(strcmp(SJ.subj, s{:})) == sum(strcmp(Mth.subj, s{:}))))
                        error('Subjects number are not matching across tasks!')
                end
                
                TAll= (index_SE) & (index_SJ) &(index_Math);
                index_subj  = strcmp(SE.subj, s{:});

                if nansum(index_subj)>0
                    prec = nansum(TAll)./nansum(index_subj);
                else
                    prec = nan;
                end % if
            end % pre_all_selfact

        end % report

        function [tout] =  numbertotalelec(obj)
            c = 0;
            for s = unique(obj.data.subj)'
                c = c+1;
                %                 tout.subj(c)  = s(:);
                tout(c)  = length(unique(obj.data.chan(categorical(obj.data.subj) == s{:})));
                %                 tout.SEselect(c) = mean(obj.table.select(categorical(obj.table.subj) == s{:} & ~obj.flag')=="Self-episodic");
                %                 tout.SJselect(c) = mean(obj.table.select(categorical(obj.table.subj) == s{:} & ~obj.flag')=="Self-judgment");
                %                 tout.view(c)        = obj.view;
                %                 tout.hemi(c)        = obj.hemi;
            end % for
        end % numbertotalelec
        function [out]    = regionals_numbers(obj, hemi, view)
            c = 0;
            sig = {};
            if strcmp(hemi, "left") && strcmp(view, "medial")
                for s = unique(obj.data.subj)'
                    c = c+1;
                    ch_c = 0;
                    for ch = unique(obj.data.chan(categorical(obj.data.subj) == s{:}))'
                        ch_c =ch_c +1;
                        sig{c}.both(ch_c)= mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:} &...
                            obj.data.Tval>0 &...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            (strcmp(obj.data.task, 'EP'))))<.05 ...
                            &...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true)&...
                            strcmp(obj.data.task, 'SJ'),:))<.05;

                        sig{c}.EP(ch_c)= mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            strcmp(obj.data.task, 'EP'),:))<.05 &...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true)&...
                            strcmp(obj.data.task, 'SJ'),:))>.05;


                        sig{c}.SJ(ch_c) =  mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true)&...
                            strcmp(obj.data.task, 'SJ'),:))<.05 &...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true)&...
                            strcmp(obj.data.task, 'EP'),:))>.05;

                    end % for ch
                end % for subj
            elseif strcmp(hemi, "right") && strcmp(view, "medial")
                for s = unique(obj.data.subj)'
                    c = c+1;
                    ch_c = 0;
                    for ch = unique(obj.data.chan(categorical(obj.data.subj) == s{:}))'
                        ch_c =ch_c +1;
                        sig{c}.both(ch_c)= mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:} &...
                            obj.data.Tval>0 &...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            (strcmp(obj.data.task, 'EP') ),:))<.05 ...
                            &...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            strcmp(obj.data.task, 'SJ'),:))<.05;

                        sig{c}.EP(ch_c)= mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            strcmp(obj.data.task, 'EP'),:))<.05&...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true)&...
                            strcmp(obj.data.task, 'SJ'),:))>.05;

                        sig{c}.SJ(ch_c) =  mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            strcmp(obj.data.task, 'SJ'),:))<.05&...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'MPFC')), obj.data.JPAnatomy, 'UniformOutput', true)&...
                            strcmp(obj.data.task, 'EP'),:))>.05;
                    end % for ch
                end % for subj
            elseif strcmp(hemi, "left") && strcmp(view, "ventral")
                for s = unique(obj.data.subj)'
                    c = c+1;
                    ch_c = 0;
                    for ch = unique(obj.data.chan(categorical(obj.data.subj) == s{:}))'
                        ch_c =ch_c +1;
                        sig{c}.both(ch_c)= mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:} &...
                            obj.data.Tval>0 &...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            (strcmp(obj.data.task, 'EP') ),:))<.05 ...
                            &...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            strcmp(obj.data.task, 'SJ'),:))<.05;

                        sig{c}.EP(ch_c)= mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            strcmp(obj.data.task, 'EP'),:))<.05 &...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true)&...
                            strcmp(obj.data.task, 'SJ'),:))>.05;

                        sig{c}.SJ(ch_c) =  mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            strcmp(obj.data.task, 'SJ'),:))<.05&...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.X < 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true)&...
                            strcmp(obj.data.task, 'EP'),:))>.05;
                    end % for ch
                end % for subj
            else
                for s = unique(obj.data.subj)'
                    c = c+1;
                    ch_c = 0;
                    for ch = unique(obj.data.chan(categorical(obj.data.subj) == s{:}))'
                        ch_c =ch_c +1;
                        sig{c}.both(ch_c)= mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:} &...
                            obj.data.Tval>0 &...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            (strcmp(obj.data.task, 'EP') ),:))<.05 ...
                            &...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            strcmp(obj.data.task, 'SJ'),:))<.05;

                        sig{c}.EP(ch_c)= mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true) &...
                            strcmp(obj.data.task, 'EP'),:))<.05 &...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true)&...
                            strcmp(obj.data.task, 'SJ'),:))>.05;

                        sig{c}.SJ(ch_c) =  mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.Tval>0 &...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true)&...
                            strcmp(obj.data.task, 'SJ'),:))<.05&...
                            mean(obj.data.Pval(categorical(obj.data.subj) == s{:} &...
                            categorical(obj.data.chan) ==  ch{:}&...
                            obj.data.X > 0 &...
                            cellfun(@(x) any(strcmp(x, 'OFC')), obj.data.JPAnatomy, 'UniformOutput', true)&...
                            strcmp(obj.data.task, 'EP'),:))>.05;
                    end % for ch
                end % for subj
            end
            out =sig;
        end % regionals_numbers
        function BHV_s = get.BHV(obj)
            % get method to load BHV structure from the json file.

            % Load JSON data from a file
            jsonStr = fileread(obj.jbpath);
            % Convert the JSON string to a struct
            BHV_s = jsondecode(jsonStr);
        end % get.BHV
        function ECoGSEEG = get.ECoGSEEG(obj)
            % get method to load structure contains number of ECoG & SEEG for each paitient.
            jsonStr = fileread(obj.jepath);
            % Convert the JSON string to a struct
            ECoGSEEG = jsondecode(jsonStr);
        end % ECoGSEEG
        function indexS = get.indexSelf(obj)
             tmp= resultiEEG.localizing(obj.data);
                indexS  = (tmp{1}.Loc_Pval<.05 & tmp{1}.Loc_Tval >0) & ...
                 ~(tmp{2}.Pval<.05 & tmp{2}.Tval <0);% not deactivated in the math   
        end

        function indexnS = get.indexNotSelf(obj)
            tmp= resultiEEG.localizing(obj.data);
            indexnS  = (tmp{1}.Loc_Pval>.05);
        end

    end % methods
end % stat_report
% $END
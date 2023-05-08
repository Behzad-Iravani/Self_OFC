classdef LMMSESJ < stat.LMM
    % CLASS NAME: LMMSESJ
    %
    % Purpose: Linear mixed-effect model for comparing self-referentially
    % activated population across SE and SJ. This
    % script is part of the codes for reproducing stats reported in
    % "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
    %
    % Properties:
    %   - coeff:        a structure containing the average of the mixed-effect model
    %                   coefficient
    %   - coeffl:       a structure containing the lower bound of the mixed-effect model
    %                   coefficient
    %   - coeffh:       a structure containing the higher bound of the  mixed-effect model
    %                   coefficient
    %   - prediction:   a vector containing the model prediction
    %   - index:        a structure containing the indices for each coefficient
    % Methods:
    %
    %   Copyright (C)  Behzad Iravani, department of neurology and neurological
    %   sciences, Stanford University. May 2023.
    %
    %   Author: Behzad Iravani
    %   behzadiravani@gmail.com
    %   Contact: behzadiravani@gmail.com
    %   Date: 05/07/2023
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        coeff   % average bootstrapped coefficient
        coeffl  % lower bound bootstrapped coefficient
        coeffh  % higher bound bootstrapped coefficient
        prediction % predicted the effect size using the model
        index  % a structure containing the indices for each coefficient
    end
    properties(Dependent)
        preprocT % preprocessed data for SE and SJ LMM
    end
    methods
        function obj = predict(obj)
            obj.prediction = struct();
            % assembeling the coeffecients:
            obj = obj.assembel_coeff();
            % predict used the bootstrapped coefficient to predict
            % fixed term
            for co = ["coeff", "coeffl", "coeffh"]
            Fixed =  obj.(co).task_EP_JPAnatomy_MPFC .* double(obj.preprocT.JPAnatomy == 'MPFC' & obj.preprocT.task == 'EP') + ...
                obj.(co).task_SJ_JPAnatomy_MPFC .* double(obj.preprocT.JPAnatomy == 'MPFC' & obj.preprocT.task == 'SJ') + ...
                obj.(co).task_EP_JPAnatomy_OFC .* double(obj.preprocT.JPAnatomy == 'OFC' & obj.preprocT.task == 'EP') + ...
                obj.(co).task_SJ_JPAnatomy_OFC .* double(obj.preprocT.JPAnatomy == 'OFC' & obj.preprocT.task == 'SJ');
            % random term
            RandEff  = zeros(height(obj.preprocT),1); % initializing the random effect term
            for s = unique(obj.preprocT.subj)' % loop over subjects
                RandEff = RandEff + obj.(co).(char(s)) .* double(obj.preprocT.subj == s);
            end % end loop over subjects
            for d = unique(obj.preprocT.Density)' % loop over electrode density
                RandEff = RandEff + obj.(co).(['d', num2str(d)]) .* double(obj.preprocT.Density == d);
            end % end loop over electrode density

            obj.prediction.(co) = Fixed + RandEff;
            end
        end % predict
        function bars(obj)
            % mdl bootstrapped model
            Coeff = cat(2,obj.mdl.Coeff)';
            %             if isempty(obj.prediction) % check if prediction is available
            %                 error('model has not been predicted yet, please run the prediction and try again!')
            %             end
            clear b ticks % clear b and ticks
            figure % open a new figure
            c = 0; % set counter
            ax1 = subplot(121); % create the axis 1
            ax2 = subplot(122); % create the axis 2
            % parse the coeff name
            obj = obj.parseCoeff();
            % create the id indices
            id = [
                find(~obj.index.OFC & obj.index.EP) % vmPFC and self-episodic
                find(~obj.index.OFC & obj.index.SJ) % vmPFC and self-judgment
                find(obj.index.OFC & obj.index.EP)  % OFC and self-episodic
                find(obj.index.OFC & obj.index.SJ)  % OFC and self-judgment
                ];
            cellfun(@(x)disp(x), obj.mdl(1).CoeffName)
            % get colors
            col = stat.LMMSESJ.colors([1,3], 2); % get two colors with two reps
            for i = 1:length(id) % loop over the predictors (task and sites)
                if i==3
                    c = 0;
                    axes(ax1)
                    title('\rmOFC', 'FontName','Arial','FontSize', 22)
                elseif i<4
                    axes(ax2)
                    title('\rmvmPFC','FontName','Arial','FontSize', 22)
                end
                c = c+1; % increment the counter
                hold on % keep the plot
                b(i) = bar(c,  quantile(Coeff(:,id(i)),.5));
                b(i).FaceColor =col(c,:);
                errorbar(c, quantile(Coeff(:,id(i)),.5), ...
                    quantile(Coeff(:,id(i)),.5)-quantile(Coeff(:,id(i)),.025), ...
                    quantile(Coeff(:,id(i)),.975)-quantile(Coeff(:,id(i)),.5), ...
                    'Color', 'k', 'Linewidth', 1, 'CapSize', 0)

                ticks(c) = median([c,c]);

            end
            ylim(ax1,[-.5,1.5])
            ylim(ax2,[-1,3])

            set(ax1, 'XTick', ticks, 'YTick', [-.5, 0,1.5], 'XTickLabel', {'SE', 'SJ', 'MTH'},...
                'TickLabelInterpreter', 'none', 'LineWidth', 2, 'FontName', 'Arial', 'FontSize', 20)
            set(ax2, 'XTick', ticks, 'YTick', [-1, 0,3], 'XTickLabel', {'SE', 'SJ', 'MTH'},...
                'TickLabelInterpreter', 'none', 'LineWidth', 2, 'FontName', 'Arial', 'FontSize', 20)

            ylabel(ax1,'Bootstrap test (95% CI)','FontName','Arial')
            ax1.TickLength(1) = .02;
            ax2.TickLength(1) = .02;
            % ax.XAxisLocation = 'origin';
            % pbaspect([1, .5, 1])

            print -dpng -r300 results\Fig1D.png
			print -dsvg results\Fig1D.svg

        end % bars

        function obj = parseCoeff(obj)
            % parseCoeff parse the coefficient's name for the given model
            obj.index.OFC = cellfun(@(x) contains(x,'OFC'), obj.mdl(1).CoeffName);
            obj.index.EP = cellfun(@(x) contains(x,'EP'), obj.mdl(1).CoeffName);
            obj.index.SJ = cellfun(@(x) contains(x,'SJ'), obj.mdl(1).CoeffName);
        end

        %----- get methods
        function pT = get.preprocT(obj)
            disp('preprocssing the input table...')
            pT = obj.preprocessing(obj.data); % preprcossing the table
            disp('done!')
        end % get method for preprocT
    end % methods
    methods(Static)
        function preprocT = preprocessing(Table)
            % prepocessing method prepares the input Table for
            % performing LMM
            % -------- find self-referential elecs
            preprocT = stat.LMMSESJ.find_self_activated_electrods(Table);

            % converting all the cell string data to categorical variable
            preprocT.subj      = categorical(preprocT.subj);
            preprocT.task      = categorical(preprocT.task);
            preprocT.JPAnatomy = categorical(preprocT.JPAnatomy);
            % remove non-self-referentially activated elecs
            preprocT(~preprocT.self_act,:) =[]; disp('removing non-self-referentially activated electrodes');
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

        function out = bootfun(Table)
            % bootfun fits mixed-effect model to Tval columns of the Table.
            % The predictors are task(2 levels: SE or SJ) and anatomical sites (i.e., JPAnatomy, 2 levels: MPFC or OFC)
            % Input:
            %       - Table:  Table containing the data including the t-values of HFB
            %                 and the labels for tasks and sites.
            % Output:
            %       - out: Structure containing the coefficients and random effect
            %       derived from LMM.
            %
            %   LMM_SEvsSJ is part of the scripts that calculates the statistics
            %   that were reported in "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
            %
            %   Copyright (C)  Behzad Iravani, department of neurology and neurological
            %   sciences, Stanford University. May 2023
            %
            %   Author: Behzad Iravani
            %   behzadiravani@gmail.com
            %   Contact: behzadiravani@gmail.com
            %   Date: 05/03/2023
            %   ------------------------------------------------------

            persistent call_
            if isempty(call_)
                call_ = 0;
            end

            call_ = call_ + 1;
            % writes to the console this number of call
            fprintf('iteration : %d\n', call_);

            % run the LMM for this iteration
            b21 = fitlme(Table, 'Tval ~ -1 + task:JPAnatomy + (1|subj) + (1|Density) ',...
                'DummyVarCoding','full'); % full dummy coding

            % store the results in out
            out.Coeff     = b21.Coefficients.Estimate;
            out.CoeffName = b21.CoefficientNames;
            out.df        = unique(b21.Coefficients.DF);
            out.rsquared  = b21.Rsquared.Ordinary;
            [~, ~, out.randomeffects_table] = randomEffects(b21);

        end % bootfun

        function Tout = find_self_activated_electrods(Tin)
            % find_self_activated_electrods finds the electrodes that are
            % activated either in EP (self-episodic/SE) or SJ
            % Input:
            %          -Tin:    a table contain the t-value comparing with math
            %                   Loc_Tval
            % Output:
            %          -Tout:   a table similar to Tin with an addition
            %                   column (self_act) that indicates if the element
            %                   was self-referentially activated
            % ---------------------------------------------------------------

            disp('finding self-referentially activated electrodes for this analysis.')
            Tout = Tin(strcmp(Tin.task, 'EP') | strcmp(Tin.task, 'SJ'),:); % initializing the output table with the input table
            warning off
            for s = unique(Tout.subj)' % loop over the unique subjects
                for ch = Tout.chan(strcmp(Tout.subj,s{:}))' % loop over the channels
                    % find self-activated electrodes for a give subject and
                    % channels
                    Tout.self_act(strcmp(Tout.subj,s{:}) & strcmp(Tout.chan, ch{:})) = ...
                        any(Tout.Loc_Tval(strcmp(Tout.subj,s{:}) & strcmp(Tout.chan, ch{:}))>0 &...
                        Tout.Loc_Pval(strcmp(Tout.subj,s{:}) & strcmp(Tout.chan, ch{:}))<=.05);
                end % inner for channel
            end % outter for subj
            warning on
        end % find_self_activated_electrods
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
            col_resource = [1.0000    0.5093    0.3465
                1.0000    0.5798    0.2505
                0.8904    0.7774         0
                0.4917    0.8684    0.2448
                0    0.8810    0.7613
                0.3902    0.8001    1.0000];
            if index> size(col_resource,1) % check if color resource has enough colors
                error('only 6 colors are available, use rep to repeate colors if more are needed!')
            end
            col = col_resource(repmat(index, rep, 1),:);
        end % colors
    end % methods Static
    % ------------------
    % private methods
    % ------------------
    methods(Access=private)
        function obj = assembel_coeff(obj)
            % assemble_coeff aggregates the coefficients of the bootstrapped
            % model.
            % Input:
            %       -mdl:    a structure containing the bootstrapped LMM
            % Output:
            %       -coeff   a structure containing average bootstrap distribution
            %       -coeffl  a structure containing low boundary of the bootstrap distribution
            %       -coeffh  a structure containing high boundary of the bootstrap distribution
            % -------------------------------------------------------------

            obj.coeffl = struct(); % initialize lower bound for coeff
            obj.coeffh = struct(); % initialize higher bound for coeff

            co = nanmean(cat(2,obj.mdl.Coeff),2); % the mean of coeff stored in co
            so = nanstd(cat(2,obj.mdl.Coeff),[],2); % the std of coeff stored in so

            obj.coeff = struct(); % initialize the coeff
            iname = 0; % set the name counter
            names = cellfun(@(x) regexprep(x, ':', '_'), obj.mdl(1).CoeffName, 'UniformOutput', false); % preprocess the field names and replace : with _ to avoid syntax problem
            for name = names % loop over the names in the mld, basically the fix terms of the model
                iname = iname +1; % increment the counter
                % storing the statistics in coeff structure
                obj.coeff.(name{:}) = co(iname);
                obj.coeffl.(name{:}) = co(iname)-1.96*so(iname);
                obj.coeffh.(name{:}) = co(iname)+1.96*so(iname);
            end % for names
            %add random
            for i = 1:length(obj.mdl) % loop over random terms
                for i2 = 1:length(obj.mdl(i).randomeffects_table.Level) % loop over levels
                    if strcmp(obj.mdl(i).randomeffects_table.Group{i2}, 'subj') % check if the random variable is subject
                        if ~isfield(obj.coeff, obj.mdl(i).randomeffects_table.Level{i2})
                            obj.coeff.(obj.mdl(i).randomeffects_table.Level{i2}) = [];
                        end
                        obj.coeff.(obj.mdl(i).randomeffects_table.Level{i2})(end+1) = obj.mdl(i).randomeffects_table.Estimate(i2);
                    else % if not subj --> density
                        namefiled = ['d', obj.mdl(i).randomeffects_table.Level{i2}]; % adding "d" to indicate this random effect for the electrodes density
                        if ~isfield(obj.coeff,namefiled)
                            obj.coeff.(namefiled) = [];
                        end
                        obj.coeff.(namefiled)(end+1) = obj.mdl(i).randomeffects_table.Estimate(i2);
                    end % if
                end % for loop over levels
            end % for loop over random terms
            % average random effect
            for i_f = fieldnames(obj.coeff)' % loop over the filed for random effecr
                if length(obj.coeff.(i_f{:}))>1
                    S = std(obj.coeff.(i_f{:}));                   % standard deviation
                    obj.coeff.(i_f{:}) = mean(obj.coeff.(i_f{:}));     % average bootstrap distribution
                    obj.coeffl.(i_f{:}) = obj.coeff.(i_f{:}) -1.96* S; % lower CI boundry
                    obj.coeffh.(i_f{:}) = obj.coeff.(i_f{:}) +1.96* S; % higher CI boundry
                end
            end % for loop over the fields for random effect

        end % assemble_coeff
    end % methods private
end % class LMMSESJ
% $END

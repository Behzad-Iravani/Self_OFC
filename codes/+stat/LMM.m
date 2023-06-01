% -*- UTF-8 -*-
classdef LMM % abstract mixed effect model object
    % CLASS NAME: LMM
    %
    % Purpose: LMM provides methods for conducting mixed-effected models
    % for assessing the HBF of OFC across task and anatomical regions. This
    % script is part of the codes for reproducing stats reported in
    % "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
    %
    % Properties:
    %   - data: a table containing the data for mixed-effect model
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
        data table   % a table containing the data
        model string % a string that defines the model to be fitted to data
        mdl % the output model of bootstrap
        dummy = 'full'
    end

    methods(Abstract)
        % Abstract method for performing preprocessing
        preprocessing(obj)
        bars(x)
        parseCoeff(obj)
    end % abstract methods

    methods
        function obj = LMM(data, model, dummy)
            % LMM Construct an instance of this abstract class
            obj.data  = data;
            obj.model = model;
            if nargin>2
            obj.dummy = dummy;
            end
        end % constructor

        function obj = bootstramp(obj,pT, np)
            % bootstrap the coefficients of LMM
            % Input:
            %       np: number of resampling
            % mdl
            %       mdl: a structure containing the summary of
            %       bootstrap statistics of LMM
            %---------------------------------------------------------
            if nargin < 3 % checks if np is provided
                np = 1e3; % default value of np is 1000
            end
            % performing bootstrapping
            disp('Bootstrapping, please wait...')
            clear obj.bootfun
% % %             s = RandStream('mlfg6331_64','Seed',1); % has substreams
% % %             opts = statset('UseParallel',true,...
% % %            'Streams',s,'UseSubstreams',true);
            opts = statset('UseParallel',true);
            obj.mdl = bootstrp(np, @obj.bootfun, pT, 'options', opts);
            disp('done!')
        end % bootstrap

        function out = bootfun(obj,Table)
            % bootfun fits mixed-effect model to Tval columns of the Table.
            % The predictors are task(2 levels: SE or SJ) and anatomical sites (i.e., JPAnatomy, 2 levels: MPFC or OFC)
            % Input:
            %       - Table:  Table containing the data including the t-values of HFB
            %                 and the labels for tasks and sites.
            % Output:
            %       - out: Structure containing the coefficients and random effect
            %       derived from LMM.
            %
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
            b21 = fitlme(Table, obj.model,...
                'DummyVarCoding',obj.dummy); % full dummy coding

            % store the results in out
            out.Coeff     = b21.Coefficients.Estimate;
            out.CoeffName = b21.CoefficientNames;
            out.df        = unique(b21.Coefficients.DF);
            out.rsquared  = b21.Rsquared.Ordinary;
            out.AIC       = b21.ModelCriterion.AIC;
            try
            [out.cf.p(1), out.cf.f(1), out.cf.df1(1), out.cf.df2(1)]     = b21.coefTest([-1, 0, 0 ,1]); % {'task_NoNeg'}    {'task_NoPos'}    {'task_YesNeg'}    {'task_YesPos'}
            out.cf.contrast(1,:) = [-1, 0, 0 ,1];
            [out.cf.p(2), out.cf.f(2), out.cf.df1(2), out.cf.df2(2)]     = b21.coefTest([0, -1, 1 ,0]);
            out.cf.contrast(2,:) = [0, -1, 1 ,0];
            [out.cf.p(3), out.cf.f(3), out.cf.df1(3), out.cf.df2(3)]     = b21.coefTest([1, -1, 1 ,-1]);
            out.cf.contrast(3,:) = [1, -1, 1 ,-1];
            [out.cf.p(4), out.cf.f(4), out.cf.df1(4), out.cf.df2(4)]     = b21.coefTest([1, 1, -1 ,-1]);
            out.cf.contrast(4,:) = [1, 1, -1 ,-1];


            end
            [~, ~, out.randomeffects_table] = randomEffects(b21);

        end % bootfun
        function report(obj, name)
            %  report is a method of object LMMSESJ that writes the bootstrap result
            %  results to .txt file.
            % -------------------------
            Coeff = cat(2,obj.mdl.Coeff)'; % get the bootstraped coefficient 
            % coefftest
            tmp = [obj.mdl.cf];
            out.f =  median(cat(1,tmp.f));
            out.p =  median(cat(1,tmp.p));
            out.df1 = ceil(median(cat(1,tmp.df1)));
            out.df2 = ceil(median(cat(1,tmp.df2)));
            out.contrast = tmp.contrast;
            out.coeffnames = obj.mdl(1).CoeffName;
            % decoding coefficient names 
            task = setdiff(fieldnames(obj.index), {'OFC'})';
            it = 0;
            for t = task
                it = it +1;
                if ismember({'OFC'},fieldnames(obj.index)')
                id(it,:) = [
                    find(obj.index.OFC & obj.index.(t{:})),find(~obj.index.OFC & obj.index.(t{:}))
                    ];
                else
                  id(it,:) = [
                    find(obj.index.(t{:}))
                    ];

                end
            end
            id = id(:);
             % open a txt file 
             disp(['writing the results to file:: results\' name '_bootResult.txt'])
            fid = fopen(['results\' name '.txt'], 'w'); % get the fild id 
            for i = 1:length(id) % loops over the predictors 
                fprintf(fid, '%s: %1.2f +/- %1.2f ci = [%1.2f, %1.2f], p = %.3g\n', obj.mdl(1).CoeffName{id(i)}, ...
                    mean(Coeff(:,id(i))), ...
                    std(Coeff(:,id(i))), ...
                    quantile(Coeff(:,id(i)),.025), quantile(Coeff(:,id(i)),.975),...
                    (1+sum(sign(mean(Coeff(:,id(i))))*Coeff(:,id(i))<=0))./(length(Coeff)+1)...
                    );
            end

            for id_ =[1,2;3,4]
                fprintf(fid, '%s-%s: %1.2f +/- %1.2f ci = [%1.2f, %1.2f], p = %.3g\n',...
                    obj.mdl(1).CoeffName{id(id_(1))},obj.mdl(1).CoeffName{id(id_(2))}, ...
                    mean(Coeff(:,id(id_(1)))-Coeff(:,id(id_(2)))), ...
                    std([Coeff(:,id(id_(1)))-Coeff(:,id(id_(2)))]), ...
                    quantile([Coeff(:,id(id_(1)))-Coeff(:,id(id_(2)))],.025), quantile([Coeff(:,id(id_(1)))-Coeff(:,id(id_(2)))],.975),...
                    (1+sum(sign(mean([Coeff(:,id(id_(1)))-Coeff(:,id(id_(2)))]))*[Coeff(:,id(id_(1)))-Coeff(:,id(id_(2)))]<=0))./(length([Coeff(:,id(id_(1)))-Coeff(:,id(id_(2)))])+1)...
                    );
            end
            fclose(fid); % close the file
            disp('done!')
            % write coeftest
            fid = fopen(['results\' name '_coefftest.json'], 'w');
            jsontxt = jsonencode(out, 'PrettyPrint',true);
            fwrite(fid, jsontxt);
            fclose(fid);

        end % report
        function n = number_trials_preConds(obj)
            % report number of dof per conds
            n = struct();
            for conds = unique(obj.preprocT.task)'
                clear tmp
                tmp = stat.average_over_subj(obj.preprocT(obj.preprocT.task ==  conds,:));
                n.(char(conds)) = sum(tmp.dof);
                fprintf('number of %s trials = %d.\n', char(conds), n.(char(conds)))
            end
        end
    end % methods
end % class LMM
% $END


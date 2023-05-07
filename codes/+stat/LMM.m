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
    end
    
    methods(Abstract) 
        % Abstract method for performing preprocessing
        preprocessing(obj)
        bootfun(x)
        bars(x)
        parseCoeff
    end % abstract methods
   
    methods
        function obj = LMM(data, model)
            % LMM Construct an instance of this abstract class 
            obj.data  = data;
            obj.model = model;  
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
           obj.mdl = bootstrp(np, @obj.bootfun, pT);
           disp('done!')
        end % bootstrap
    end % methods
end % class LMM
% $END


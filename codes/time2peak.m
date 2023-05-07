%-*- UTF8 -*-
classdef time2peak
   % CLASS NAME: time2peak
    %
    % Purpose: time2peak provides methods for reporting the individual
    % onset latency for HFB self-referntial response
    %
    % "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
    %
    % Properties:
    %   - SE: a table containing the summary of data for self-episodic
    %   - SJ: a table containing the summary of data for self-judgment
    % Methods:
    %
    % Author: Behzad Iravani
    % Contact: behzadiravani@gmail.com
    % Date: 05/06/2023
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        SE table % self-episodic data
        SJ table % self-judgment data
    end

    methods
        function obj = time2peak(SE, SJ)
            % time2peak Construct an instance of this class
            obj.SE = SE;
            obj.SJ = SJ;
        end % constructor 

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.SE + inputArg;
        end
    end
end
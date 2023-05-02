classdef resultiEEG < stat_report
    % CLASS NAME: resutiEEG
    %
    % Purpose: resutiEEG provides methods for recreating the results and figures reported in 
    % "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX ",
    % 
    %
    % Properties:
    %  subclass of stat_report    
    %  
    %
    % Methods:
    %
    % Author: Behzad Iravani
    % Contact: behzadiravani@gmail.com
    % Date: 05/02/2023
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties(Dependent)
        surf_ % surf object for plotting brain and electrodes 
    end
    
methods
    function obj = resultiEEG(data, jbpath, jepath)
       % constructor method for creating instance of resultiEEG
        obj@stat_report(data, jbpath, jepath)
    end

end % methods

end
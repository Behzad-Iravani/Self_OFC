% -*- UTF 8 -*-
function [s, D, flag_remove, label] = plot_electrodes_projected_to_brain(electrodeMNI, MNI_surface, hemi, view, JPlable, szx)
% plot_electrodes_projected_to_brain plots the electrodes on the brain
% the electrodes selected to be plotted depends on the view 
% Input:
%       - electrodeMNI: MNI coordinates of the electrodes
%       - MNIsurface:   surface object 
%       - hemi:         determines which side to be plotted (i.e., left or right) 
%       - view:         view to be plotted see bellow for more information 
%       - JPlabel:      the anatomical lable of each electrodes 
%       - szx:          the size of electrodes, scalar
% Output:
%       - s:            scatter plot object for each electrodes 
%       - D:            depth of the electrode from the surface of the
%                       cortex.
%       - flag_remove:  if the electrode should not be plotted for the
%                       given view
%       - label:        the anatmoical lable of the electrode.
%
%   plot_electAcitivity is part of the scripts that recreates the plots and
%   that was reported in "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
%
%   Copyright (C)  Behzad Iravani, department of neurology and neurological
%   sciences, Stanford University. May 2023
%
%   Author: Behzad Iravani
%   behzadiravani@gmail.com
%   Contact: behzadiravani@gmail.com
%   Date: 05/02/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 0;
switch view % which view to plot the brain 
    case 'medial'
        projection = 2;  % loc_view(90, 0)
    case 'lateral'
        projection = 2;
    case 'anterior'
        projection = 1;
    case 'posterior'
        loc_view(0,0)
    case 'ventral'
        projection = 2;
    case 'temporal'
        projection = 2;
    case 'dorsal'
        projection = 2;
    case 'latero-ventral'
        projection = 2;
    case 'medio-dorsal'
        projection = 2;
    case 'medio-ventral'
        projection = 2;
    case 'medio-posterior'
        projection = 2;
    case 'medio-anterior'
        projection = 2;
    case 'frontal'
        projection = 1;
    case 'parietal'
        projection = 2;
end


for i  = 1:size(electrodeMNI,1)
    [dist,idx] = min(sum((electrodeMNI(i,:)-MNI_surface.(hemi{:}).cortex.vert).^2,2));

    if projection == 1
        x = electrodeMNI(i,1);
        y = max(MNI_surface.(hemi{:}).cortex.vert(:,2));
        z = electrodeMNI(i,3);
    elseif projection == 2
        switch view
            case 'medial'
                x = 0;
                y = electrodeMNI(i,2);
                z = electrodeMNI(i,3);
                if ~isempty(JPlable{i})
                    if iscell(JPlable{i})
                    flag_remove(i) = ~any(cellfun(@(x) strcmp(x,'MPFC') | strcmp(x,'MFC'), JPlable{i}));
                    else
                        flag_remove(i) = ~(strcmp(JPlable{i},'MPFC') | strcmp(JPlable{i},'MFC'));
                    end
                else
                    flag_remove(i) = false;
                end
                if ~flag_remove(i)

                    label{i} = 'MFC';
                else
                    label{i} = nan;
                end
            case 'ventral'
                x = electrodeMNI(i,1);
                y = electrodeMNI(i,2);
                z = min(MNI_surface.(hemi{:}).cortex.vert(:,3));
                if ~isempty(JPlable{i})
                    if iscell(JPlable{i})
                        flag_remove(i) = ~any(cellfun(@(x) strcmp(x,'OFC'), JPlable{i}));
                    else
                       flag_remove(i) = ~strcmp(JPlable{i},'OFC');
                    end
                else
                    flag_remove(i) = false;
                end
                if ~flag_remove(i)

                    label{i} = 'OFC';
                else
                    label{i} = nan;
                end
        end

    else
        x = electrodeMNI(i,1);
        y = electrodeMNI(i,2);
        z = min(MNI_surface.(hemi{:}).cortex.vert(idx,3));
    end

    if strcmp(hemi{:},'right') && electrodeMNI(i,1)>0 
        c = c+1;
        s(c) = scatter3(x,y,z,szx,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.25);
        D(c) = dist;
    elseif strcmp(hemi{:},'left') && electrodeMNI(i,1)<0
        c = c+1;
        s(c) = scatter3(x,y,z,szx,'filled','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.25);
        D(c) = dist;
   
    end

end
% $ END

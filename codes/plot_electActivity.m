%-*- UTF8 -*-
function plot_electActivity(MNI_surface, hemi, right, left, col)
% plot_electActivity plots color coded the electrodes on the brain surface
% Input:
%       - MNI_surface: surf_ object contains the surface and electrode data.
%       - hemi:        cell that defines which hemisphere to plot. e.g., {'right'}
%       - right:       
%       - left
%       - col
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


% load MNI surafce of hemi
[vertex_coords, faces] = MNI_surface.read_surf(MNI_surface.(hemi{:}).path);

if ~isempty(find(faces == 0)) % checks if the indexing starts with 0
    warning('indxing statrs at zero adding 1 to faces')
    faces = faces +1;
end
% assing the vetrxes and face to surface object 
MNI_surface.(hemi{:}).cortex.vert  = vertex_coords;
MNI_surface.(hemi{:}).cortex.tri   = faces;
% computes density map and correcting
MNI_surface = MNI_surface.densityMAP(hemi{:}); % computes the density of electrodes

% plotting  brain
for views = {'ventral', 'medial'} % loop over the views 
    figure % open a new figure per each view & hemispheres 
    subplot(1,10,[1:9]);
    MNI_surface.plot_brain(hemi, [.95 .95 .85]); % plot the brain 
    alpha(1)
    hold on
    if strcmp(hemi{:}, 'right') % checks the hemi side 
        if any(startsWith(strtrim(right.Properties.VariableNames), 'select'))
            actvie_elec = right.select ;
        else % not categorical
            actvie_elec = right.Loc_Tval;
        end % if active 
       
        % project electrodes on the brain 
        [s, D, flag, ~] = misc.plot_electrodes_projected_to_brain(...
            [right.X, right.Y, right.Z], MNI_surface, hemi, views{:}, right.JPAnatomy, 45);
    else % left side ----
        if any(startsWith(strtrim(left.Properties.VariableNames), 'select'))
            actvie_elec = left.select ;
        else % not categorical 
            actvie_elec = left.Loc_Tval;
        end
       
        [s, D, flag, ~] = misc.plot_electrodes_projected_to_brain(...
            [left.X, left.Y, left.Z], MNI_surface, hemi, views{:}, left.JPAnatomy, 45);

    end
    D = .3*(D/max(D))+.3; % rescaling the depth parameter for visual purposes 
    misc.reorgnize_electrodes(s, D, actvie_elec, flag, col, right, left, hemi, views, false)
    misc.plot_setting(hemi{:},views{:}) % adjust the plot setting given the view 

    axis equal
    axis tight
    axis off
    camlight
   %  add colorbar axes
    axbar  = subplot(1,10,10); % get the axes for the colorbar 
    patch([0, 1, 1, 0], [0, 0, 1, 1], [0 0 1 1])
    colm = [flipud(arrayfun(@(x) misc.adjust_sat(col(1), x), 0:.05:.5, 'UniformOutput', false)'),...
        arrayfun(@(x) misc.adjust_sat(col(2), x), 0:.05:.5, 'UniformOutput', false)'];

    colormap(axbar, flipud(cat(1, colm{:})))
    set(gca, 'Xtick', [] , 'Ytick', [0, .5 1], 'YTickLabel', [-round(abs(max([right.Loc_Tval; left.Loc_Tval]))),...
        0,...
        round(max(abs([right.Loc_Tval; left.Loc_Tval])))], 'FontName', 'Arial', 'FontSize', 12)
    % ----------------------
    print(['results\Fig1b_', hemi{:}, '_' , views{:}, '.png'], '-r600', '-dpng') % print the plot to file 
end % for view
end 
    %-----
    % $ END

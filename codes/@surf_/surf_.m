classdef surf_
    % CLASS NAME: surf_
    %
    % Purpose: surf_  represents brain surface data and electrode activity.
    %
    % Properties:
    %   - left:  store information about the left hemisphere of the brain
    %   surface.
    %   - right: store information about the right hemisphere of the brain surface.  
    %   - Elecpos: matrix that stores electrode positionsMatrix that
    %   stores electrode positions.
    %   - ElectActivity: matrix that stores electrode activity.   
    %     
    % Methods:
    %   - report(what): report the info based on the input (what),
    %   and behvaioral .json file see below for more info.
    %
    % Author: Behzad Iravani
    % Contact: behzadiravani@gmail.com
    % Date: 05/02/2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        left            = struct('path',  '@surf_\dat\lh.surf',...
                                 'cortex',[], ...
                                'CortexElectrodeDensity', [],...
                                'ElectrodDensity' , [],...
                                'ElecActivity',[],...
                                'ElecPOS', [],...
                                'Corrected_Activity', [],...
                                'spatialFilt', []);

        right           = struct('path',  '@surf_\dat\rh.surf',...
                                'cortex',[], ...
                                'CortexElectrodeDensity', [],...
                                'ElectrodDensity' , [],...
                                'ElecActivity',[],...
                                'ElecPOS', [], ...
                                'Corrected_Activity', [],...
                                'spatialFilt', []);
        ElectPos        = [];
        ElectActivity   = [];

    end % properties 

    methods
        function obj =  densityMAP(obj,hemi)
            % densityMAP calculates the density map of electrode activity on the brain surface
            % Input:
            %        - hemi:  cell contains the hemisphere lable to plot (e.g., {'right'})
            % -----------------------------------------------------------------------------

            if  ~isempty(obj.ElectPos) % Check if electrode positions are available
                % Define a 3D Gaussian kernel to calculate the density map
                FWHM  = 6; % Full width at half maximum of the Gaussian kernel, in mm
                sigma = FWHM./2.355; % standard deviation 
                G = @(x,y,z) (1/(2*pi*sigma.^3))*exp(-(x.^2+y.^2+z.^2)./(2*sigma.^3));
                % Select electrodes and their activity for the specified hemisphere
                switch hemi
                    case 'right'
                        obj.(hemi).ElecPOS = obj.ElectPos(obj.ElectPos(:,1)>0,:);
                        obj.(hemi).ElecActivity = obj.ElectActivity(obj.ElectPos(:,1)>0,:);

                    case 'left'
                        obj.(hemi).ElecPOS = obj.ElectPos(obj.ElectPos(:,1)<0,:);
                        obj.(hemi).ElecActivity = obj.ElectActivity(obj.ElectPos(:,1)<0,:);
                end % switch 
                % Calculate the density map for each vertex on the brain surface
                parfor i=1:size(obj.(hemi).ElecPOS,1)
                    fprintf('observation %d of %d\n', i, size(obj.(hemi).ElecPOS,1))
                    prob(i, :) = G(obj.(hemi).ElecPOS(i,1) - obj.(hemi).cortex.vert(:,1), ...
                        obj.(hemi).ElecPOS(i,2) - obj.(hemi).cortex.vert(:,2), ...
                        obj.(hemi).ElecPOS(i,3) - obj.(hemi).cortex.vert(:,3))./...
                        sum(G(obj.(hemi).ElecPOS(i,1) - obj.(hemi).cortex.vert(:,1), ...
                        obj.(hemi).ElecPOS(i,2) - obj.(hemi).cortex.vert(:,2), ...
                        obj.(hemi).ElecPOS(i,3) - obj.(hemi).cortex.vert(:,3)));
                end
                % Aggregate the probabilities for each vertex
                prob(isnan(prob)) = 0;
                prob_agg = sum(prob);
                % Normalize the density map and store it in the object
                obj.(hemi).CortexElectrodeDensity= 1e2*prob_agg./sum(prob_agg);
                % Map the cortex density values to the electrodes and calculate the corrected activity
                if length(prob_agg)>2
                    obj = obj.map_cortex_dens_2_elect(hemi);
                    % corrected activity for density
                    [~, ~, obj.(hemi).Corrected_Activity] = regress(obj.(hemi).ElecActivity, ...
                        [obj.(hemi).ElectrodDensity]);
                else
                    warning('no density correction due to low number of electrodes.')
                    obj.(hemi).Corrected_Activity = obj.(hemi).ElecActivity;
                end
            else
                warning('No electrodes are available, please provide ElectrodPos in MNI space')
            end
        end % densityMAP

        function obj =  WeightedMAP(obj,hemi)
             % WeightedMAP calculates the weighted heatmap of electrode activity on the brain surface
             % Input:
             %        - hemi:  cell contains the hemisphere lable to plot (e.g., {'right'})
             % -----------------------------------------------------------------------------

            if  ~isempty(obj.ElectPos)
                % Define a 3D Gaussian kernel to calculate the weighted map
                FWHM  = 6; % Full width at half maximum of the Gaussian kernel, in mm
                sigma = FWHM./2.355;% standard deviation
                G = @(x,y,z) exp(-(x.^2+y.^2+z.^2)./(2*sigma.^3));

                parfor i=1:size(obj.(hemi).ElecPOS,1)
                    fprintf('observation %d of %d\n', i, size(obj.(hemi).ElecPOS,1))
                    act(i, :) = obj.(hemi).Corrected_Activity(i) *G(obj.(hemi).ElecPOS(i,1) - obj.(hemi).cortex.vert(:,1), ...
                        obj.(hemi).ElecPOS(i,2) - obj.(hemi).cortex.vert(:,2), ...
                        obj.(hemi).ElecPOS(i,3) - obj.(hemi).cortex.vert(:,3));
                end
                 % Aggregate the probabilities for each vertex
                act(isnan(act)) = 0;
                obj.(hemi).spatialFilt = act;
                act_agg = sum(act);
                 % Normalize the density map and store it in the object
                obj.(hemi).CortexElectrodeActivity= act_agg;

            else
                warning('No electrodes are available, please provide ElectrodPos in MNI space')
            end
        end % WeightedMAP

        function obj = map_cortex_dens_2_elect(obj, hemi)
            % map_cortex_dens_2_elect maps the cortex densisty to each
            % contact. 
            % Input:
            %        - hemi:  cell contains the hemisphere lable to plot (e.g., {'right'})
            % ---------------------------------------------------------------------------

            [~,idx]= min(... Find the index of the closest vertex in the cortical surface for each electrode.
                (obj.(hemi).cortex.vert(:,1)-obj.(hemi).ElecPOS(:,1)').^2 + ...
                (obj.(hemi).cortex.vert(:,2)-obj.(hemi).ElecPOS(:,2)').^2 + ...
                (obj.(hemi).cortex.vert(:,3)-obj.(hemi).ElecPOS(:,3)').^2  ...
                ,[],1); % Compute the Euclidean distance between each vertex in the cortical surface and each electrode.
            % Assign the cortical density value of the closest vertex to each electrode, and save it in the ElectrodDensity
            obj.(hemi).ElectrodDensity =  obj.(hemi).CortexElectrodeDensity(idx)'; 
        end % map_cortex_dens_2_elect

        function obj = map_ctx_act_2_elect(obj, hemi)
             % map_cortex_act_2_elect maps the cortex activity to each
             % contact.
             % Input:
             %        - hemi:  cell contains the hemisphere lable to plot (e.g., {'right'})
             % ---------------------------------------------------------------------------

             % Define a 3D Gaussian kernel to calculate the weighted map
             FWHM  = 6; % Full width at half maximum of the Gaussian kernel, in mm
             sigma = FWHM./2.355;% standard deviation
             
             G = @(r) exp(-(r.^2)./(2*sigma.^3));
            for i= 1:length(obj.(hemi).ElecPOS)
                [rm,idx]= min(...
                    (obj.(hemi).cortex.vert(:,1)-obj.(hemi).ElecPOS(i,1)').^2 + ...
                    (obj.(hemi).cortex.vert(:,2)-obj.(hemi).ElecPOS(i,2)').^2 + ...
                    (obj.(hemi).cortex.vert(:,3)-obj.(hemi).ElecPOS(i,3)').^2  ...
                    ,[],1);

                obj.(hemi).ElecActivity(i, :)=[...
                    obj.(hemi).CortexElectrodeActivity(idx).*G(rm)];
            end % for
        end % map_ctx_act_2_elect

        function [ax, trs]  = viz_electrode_density_map(obj,hemi, col)
            % viz_electrode_density_map visualizes the electrode density map on the cortex surface
            % Input:
            %        - hemi:  cell contains the hemisphere lable to plot (e.g., {'right'})
            %        - col:  colormap matrix Nx3
            % -------------------------------------------------------------
            
            clf % clears the current figure
            ax = axes(); % creates a new set of axes in the current figure and returns the handle 
            trs = trisurf(obj.(hemi{:}).cortex.tri,...
                obj.(hemi{:}).cortex.vert(:,1),obj.(hemi{:}).cortex.vert(:,2),obj.(hemi{:}).cortex.vert(:,3),...
                obj.(hemi{:}).CortexElectrodeDensity); % triangulation of the brain surface.
            trs.EdgeColor = 'none'; % sets the edge color of the surface to none.
            trs.FaceAlpha = .75; % sets the transparency of the surface.

            colormap(ax, col) % sets the colormap of the axes
        end

        function [ax, trs]  = plot_brain(obj,hemi, col)
            % plot_brain visualizes the cortex surface
            % Input:
            %        - hemi:  cell contains the hemisphere lable to plot (e.g., {'right'})
            %        - col:  colormap matrix Nx3
            % -------------------------------------------------------------
            trs = trisurf(obj.(hemi{:}).cortex.tri,...
                obj.(hemi{:}).cortex.vert(:,1),obj.(hemi{:}).cortex.vert(:,2),obj.(hemi{:}).cortex.vert(:,3));% triangulation of the brain surface.
            trs.EdgeColor = 'none'; % sets the edge color of the surface to none.
            trs.FaceAlpha = .75;  % sets the transparency of the surface.

            colormap(gca, col) % sets the colormap of the axes
            ax = get(gca);

        end

        function [ax, trs]  = viz_electrode_act_map(obj,hemi, col)
            % viz_electrode_act_map visualizes the electrode activity map on the cortex surface
            % Input:
            %        - hemi:  cell contains the hemisphere lable to plot (e.g., {'right'})
            %        - col:  colormap matrix Nx3
            % -------------------------------------------------------------
            
            clf % clears the current figure
            ax = axes(); % creates a new set of axes in the current figure and returns the handle
            trs = trisurf(obj.(hemi{:}).cortex.tri,...
                obj.(hemi{:}).cortex.vert(:,1),obj.(hemi{:}).cortex.vert(:,2),obj.(hemi{:}).cortex.vert(:,3),...
                obj.(hemi{:}).CortexElectrodeActivity);% triangulation of the brain surface.
            trs.EdgeColor = 'none'; % sets the edge color of the surface to none.
            trs.FaceAlpha = .75;  % sets the transparency of the surface.

            colormap(gca, col) % sets the colormap of the axes
        end
    end
    methods (Static)
        % method to read bytes
        byte = fd3(fid)
        % method to read the .surf file
        [vertex_coords, faces]  = read_surf(path)
    end % Static methods

end % surf_ calss
% $ END






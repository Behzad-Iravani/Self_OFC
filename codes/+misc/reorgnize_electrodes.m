% -*- UTF-8 -*-
function reorgnize_electrodes(s, D, actvie_elec, flag, col, right, left, hemi, views, changeSEEG_C)
% reorgnize_electrode makes sure that active electrodes appear on the front
%-------------------------------------------------------------------------
if strcmp(hemi{:}, 'right')
    if ismember({'Tn'}, right.Properties.VariableNames)
        p = right.Tn;
    else
        p = right.Tval;
    end
else
    if ismember({'Tn'}, left.Properties.VariableNames)
        p = left.Tn;
    else
        p = left.Tval;
    end
end

% ----------------------ADJUSTING ELECTRODES ---------------------------
rm = []; %
if iscategorical(actvie_elec) % checks if the coloring is based on the categories
    for i = 1 :length(s)
      
        if ~isnan(D(i)) && ~isempty(s(i))
            if ~flag(i) &&(actvie_elec(i) == "Self")
                if ~changeSEEG_C
                    s(i).MarkerFaceColor =  misc.adjust_sat(col(1), p(i));%col(1,:);
                end
                s(i).MarkerFaceAlpha = .85;
            elseif ~flag(i) && (actvie_elec(i) == "Math")
                if ~changeSEEG_C
                    s(i).MarkerFaceColor = misc.adjust_sat(col(2), -p(i));%col(1,:);
                end
                s(i).MarkerFaceAlpha = .85;

            elseif ~flag(i) && actvie_elec(i)=="both"
                if ~changeSEEG_C
                    s(i).MarkerFaceColor = misc.adjust_sat([255,29,206]/255, p(i));%col(1,:);
                end
                s(i).MarkerFaceAlpha = .85;

            elseif ~flag(i) && (actvie_elec(i) == "none")
                if ~changeSEEG_C
                    s(i).MarkerFaceColor = [.9 .9 .9];
                end
                s(i).MarkerFaceAlpha = .65;
            else
                delete(s(i))
                rm = [rm,i];
            end
        end
    end
else % not categorical ---------
    for i = 1:length(s)
        %         s(i).SizeData = 20;
        if ~isnan(D(i)) && ~isempty(s(i))
            if ~flag(i) && actvie_elec(i)>0
                if ~changeSEEG_C
                    s(i).MarkerFaceColor =  misc.adjust_sat(col(1,:), p(i));%col(1,:);
                end
                s(i).MarkerFaceAlpha = .85;
            elseif ~flag(i) && actvie_elec(i)<0
                if ~changeSEEG_C
                    s(i).MarkerFaceColor = misc.adjust_sat(col(2,:), -p(i));%col(1,:);
                end
                s(i).MarkerFaceAlpha = .85;
            else
                delete(s(i))
                rm = [rm,i];
            end
        end
    end
end
% retriving the electrode types
% get the kept elect type
kept = eval(['setdiff(1:height(', hemi{:}, '),[' num2str(rm),'])']);
ElecType = eval(['{', hemi{:}, '.elecType{kept}}']);

s(rm)           = [];
flag(rm)        = [];
actvie_elec(rm) = [];
% loop over the electrodes and adjust the coordinate for visual
% purposes
for i = 1:length(s) % loop over the electrodes and bring the active one to front for viusal purposes
    if strcmp(views{:},'ventral') % checks the view side
        if iscategorical(actvie_elec) % brings the active to front to be visible
            s(i).ZData = s(i).ZData * 1.25 - double(strcmp(ElecType{i},'SEEG'))*20 ... % braing SEEG to front
                - double(strcmp(ElecType{i},'SEEG')|  (~flag(i) & ~(actvie_elec(i) == "none")))*10;
        else
            s(i).ZData = s(i).ZData * 1.25 - double(strcmp(ElecType{i},'SEEG'))*20 ... % braing SEEG to front
                - double(strcmp(ElecType{i},'SEEG') | (~flag(i) & abs(actvie_elec(i))>1))*10; % bring the active ECoG front but berfore SEEG
        end
    elseif strcmp(views{:},'medial')
        s(i).SizeData = 85; % increase the size of electrodes for medial view to be visually nice
        if strcmp(hemi{:}, 'right')  % checks the side
            if iscategorical(actvie_elec) % brings the active to front to be visible
                s(i).XData = s(i).XData - 10 -  (~flag(i) & ~(actvie_elec(i) == "none"))*10;
            else % not active
                s(i).XData = s(i).ZData * 1.25 -10 - strcmp(ElecType{i},'SEEG')* 12 ... % braing SEEG to front
                    -(strcmp(ElecType{i},'SEEG')| (~flag(i) & abs(actvie_elec(i))>1))*10; % bring the active ECoG front but berfore SEEG
            end % if active
        end
    else % neither ventral nor medial
        if iscategorical(actvie_elec) % brings the active to front to be visible
            s(i).XData = s(i).XData + 10 + (~flag(i) & ~(actvie_elec(i) == "none"))*10;
        else % not active
            s(i).XData = s(i).ZData * 1.25 +10  + strcmp(ElecType{i},'SEEG')* 12 ... % braing SEEG to front
                +(strcmp(ElecType{i},'SEEG')| (~flag(i) & abs(actvie_elec(i))>1))*10; % bring the active ECoG front but berfore SEEG
        end % if active

    end % if view
    if strcmp(ElecType{i},'SEEG') % seeg change marker to square
        if changeSEEG_C
            s(i).Marker     = 'square'; 
            s(i).SizeData = s(i).SizeData*1.50;
            s(i).MarkerEdgeColor = [.25,.25,.25];
        end
    end
end % if electrodes as scatter plot (s)
end % reorgnize_electrodes
% $ END
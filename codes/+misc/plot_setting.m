function plot_setting(hemi,viewside)

lighting gouraud; %play with lighting...

material([.2 1 .2 10 1]);

% hemi- 'l' for left hemipshere, 'r' for right hemisphere
% views are labelled. A dash between views (Ex: medio-ventral) means that the view is in between the medial and ventral views

if strcmp(hemi,'left')
    switch viewside
        case 'medial'
            loc_view(90, 0)
        case 'lateral'
            loc_view(270, 25)
        case 'anterior'
            loc_view(180,0)
        case 'posterior'
            loc_view(0,0)
        case 'ventral' 
            loc_view(180,270)
        case 'temporal'
            loc_view(270,-80)
        case 'dorsal'
            loc_view(0,90)
        case 'latero-ventral'
            loc_view(270,-45)
        case 'medio-dorsal'
            loc_view(90,45)
        case 'medio-ventral'
            loc_view(90,-45)
        case 'medio-posterior'
            loc_view(45,0)
        case 'medio-anterior'
            loc_view(135,0)
        case 'frontal'
            loc_view(-120,10)
        case 'parietal'
            loc_view(-70,10)
    end
%     set(l,'Position',[-1 0 1])   
elseif strcmp(hemi,'right')
    switch viewside
        case 'medial'
            loc_view(270, 0)
        case 'lateral'
            loc_view(90, 25)
         case 'anterior'
            loc_view(180,0)
        case 'posterior'
            loc_view(0,0)
        case 'dorsal'
            loc_view(0,90)
        case 'ventral'
            loc_view(180,270)
        case 'temporal'
            loc_view(90,-80)
       case 'latero-ventral'
            loc_view(90,-45)
        case 'medio-dorsal'
            loc_view(270,45)
        case 'medio-ventral'
            loc_view(270,-45)
        case 'medio-posterior'
            loc_view(315,0)
        case 'medio-anterior'
            loc_view(225,0)  
        case 'frontal'
            loc_view(120,10)
        case 'parietal'
            loc_view(70,10)
    end
    
else
    error('hemisphere should be either leftl or right')
%     set(l,'Position',[1 0 1])     
end



set(gcf,'Renderer', 'zbuffer')
set(gcf,'color','w')
end 
% $ END
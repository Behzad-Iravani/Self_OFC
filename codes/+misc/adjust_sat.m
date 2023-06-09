% -*- UTF-8 -*-
function col = adjust_sat(col, p)
% adjust_sat adjusts the colors staturation based on the parameter p (p-value)
% Input:
%       - col: hexcolor values 
%       - p: parameter to adjust the saturation for example p-value 
%
%   adjust_sat is part of the scripts that recreates the plots 
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

col = misc.hex2rgb(col); % converts hex color to rgb values
hsv = rgb2hsv(col); % converts the rgb value to hsv colorspace to access the saturation values 

p(p<0) = 0; % no negtaive coloring // below baseline

hsv(2) = 2*((1./(1+exp(-3*abs(p))))-.5);% a sigmoid for smooth saturation transition for visual purposes. 
col = hsv2rgb(hsv); % converts back to rgb colorspace 

end
% $ END
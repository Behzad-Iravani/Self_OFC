% -*- UTF-8 -*-
function errorbar_bootstrap(Coeff,id, id_, c, H)
% errorbar_bootstrap calculates and plots the error bar for the input Coeff
% Input:
%       Coeff:       a NxM matrix of bootstrapped coefficients, N, numeber of resampling and M number of predictors  
%       id:          a vecotor contaning the order of predictors to plot
%       id_:         a scalar index indicating the current error bar to be
%                    plotted
%         c:         a scalar inidcating the position of current error bar to be plotted
%         H:         a string that detemines the layout of the bar plot
% Author: Behzad Iravani
% Email: behzadiravani@gmail.com        
% May 9th, 2023 Philladelphia, PA
%--------------------------------------------------------------------
m = quantile(Coeff(:,id(id_,:)),.5); % median
lowerbound = abs(quantile(Coeff(:,id(id_,:)),.025)-m); 
upperbound = abs(quantile(Coeff(:,id(id_,:)),.975)-m); 
if strcmp(H,'off') % check the layout of the plot
errorbar(c, m,...
    lowerbound, upperbound,...
    'Color', 'k', 'Linewidth', 1, 'CapSize', 0)
else
errorbar(m, c,...
    lowerbound, upperbound, 'horizontal',...
    'Color', 'k', 'Linewidth', 1, 'CapSize', 0)
end
end % errorbar_bootstrap
% $END
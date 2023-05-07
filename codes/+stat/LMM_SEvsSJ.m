% -*- UFT-8 -*-
function myout = LMM_SEvsSJ(Table)
% LMM_SEvsSJ fits mixed-effect model to Tval columns of the Table. 
% The predictors are task(2 levels: SE or SJ) and anatomical sites (i.e., JPAnatomy, 2 levels: MPFC or OFC)
% Input:
%       - Table:  Table contianig the data including the t-values of HFB
%                 and the lables for tasks and sites.
% Output:
%       - myout: Structure containing the coeffecients and randome effect
%       derived from LMM.
%   
%   LMM_SEvsSJ is part of the scripts that calculates the statistics
%   that were reported in "SELF-REFERENTIAL PROCESSING IN NEURONAL POPULATIONS OF VENTROMEDIAL AND ORBITOFRONTAL CORTEX "
%   
%   Copyright (C)  Behzad Iravani, department of neurology and neurological
%   sciences, Stanford University. May 2023
%
%   Author: Behzad Iravani
%   behzadiravani@gmail.com
%   Contact: behzadiravani@gmail.com
%   Date: 05/03/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent call_
if isempty(call_)
    call_ = 0;
end

call_ = call_ + 1;
% writes to the consoule this number of call 
fprintf('iteration : %d\n', call_);

% run the LMM for this iteration 
b21 = fitlme(Table, 'Tval ~ -1 + task:JPAnatomy + (1|subj) + (1|Density) ',...
    'DummyVarCoding','full');%,Distribution='Binomial');

% store the results in myout 
myout.Coeff     = b21.Coefficients.Estimate;
myout.CoeffName = b21.CoefficientNames;
myout.df        = unique(b21.Coefficients.DF);
myout.rsquared  = b21.Rsquared.Ordinary;
[~, ~, myout.randomeffects_table] = randomEffects(b21);

end
% $END
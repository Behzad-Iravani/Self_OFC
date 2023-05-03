% -*- UFT-8 -*-
function [t, all, n, s, ci] = t_varaible_length(T, lock_type,s)
% t_varaible_length calculates the statistics for table T

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = stat.average_over_task(T); % average over the task

mxt = max(cellfun(@length,T.avg));
all = zeros(height(T), mxt);
tmp = zeros(1, mxt);
n   = zeros(1, mxt);

baseline =[];
if strcmp(lock_type, 'resp')
for h=1:height(T)
    tmpa = zeros(1, mxt);
    tmpa(end-length(T.avg{h})+1:end) = T.avg{h};
    baseline(h,:) = mean(T.avg{h}(T.time{h}>0)); % the baseline for the response lock is after zero 
    all(h,:) = movmean(tmpa ,floor(s*1e3));
    tmp = tmp + tmpa ;
    
    tmpn = zeros(1, mxt);
    tmpn(end-length(T.avg{h})+1:end) = 1;
    n = n + tmpn;
end
elseif strcmp(lock_type, 'stim')
for h=1:height(T)
    tmpa = zeros(1, mxt);
    tmpa(1:length(T.avg{h})) = T.avg{h};
    baseline(h,:) = mean(T.avg{h}(T.time{h}<0));
    all(h,:) = movmean(tmpa ,floor(s*1e3));
    tmp = tmp + tmpa;

   
    tmpn = zeros(1, mxt);
    tmpn(1:length(T.avg{h})) = 1;
    n = n + tmpn;
end

end
m = tmp./n;
s  = sum((all-repmat(m, size(all,1),1)).^2)./(n-1);

% % % sb = (tmpb-mean(tmpb)).^2./(n-1);
% tb = (tmpb -mean(tmpb))./sqrt(sb);
t = sqrt(n).*(m-mean(baseline))./sqrt(s);
t_value = tinv(1 - .05/2, n-1);
ci(1,:) = (m -mean(baseline)- t_value.*sqrt(s./n));
ci(2,:) = (m -mean(baseline)+ t_value.*sqrt(s./n));






function [mean_lfp,std_lfp,num_counted] = combineVariableLengthLfps(lfps)
if iscell(lfps)
  lfp_mat = cellArray2mat(lfps,'single');
else
  lfp_mat = lfps;
end
num_counted = sum(~isnan(lfp_mat),2);
idxs = num_counted > 1;

mean_lfp = mean(lfp_mat,2,'omitnan');
mean_lfp = mean_lfp(idxs);
std_lfp = std(lfp_mat,0,2,'omitnan');
std_lfp = std_lfp(idxs) ./ (sqrt(num_counted(idxs))-1);
end
function [binned_cohs, num_in_bins, mean_cohs, std_err] = groupCohByDist(cohs, dist_edges)
if ndims(cohs) == 3
  tmp = zeros([size(cohs),1]); % because reshape(cohs, ..., 1) doesn't work
  tmp(:,:,:,1) = cohs;
  cohs = tmp;
  clear tmp;
end
s = size(cohs);
assert(s(2) == 96 && s(3) == 96);
num_trials = s(end);
% collapse inner dims for easy indexing
cohs = reshape(cohs, [s(1), 96*96, num_trials]);

idxs = getCohDistIdxs(dist_edges);
nbins = numel(idxs);

mean_cohs = zeros(s(1),num_trials,nbins);
num_in_bins = zeros(1,nbins);
binned_cohs = cell(1,nbins);
std_err = zeros(size(mean_cohs));
for i=1:nbins
  num_in_bins(i) = size(idxs{i},1);
  % make 96x96 matrix where pos i,j is 1 iff elecs i and j should be used
  idx_mat = logical(accumarray(idxs{i},1,[96,96]));
  binned_cohs{i} = cohs(:,idx_mat(:),:);
end
if nargout > 2 % this can be expensive, so there's a way to avoid this calculation
  for ii=1:nbins
    [mean_cohs(:,:,ii), std_err(:,:,ii)] = calcMeanAndStd(binned_cohs{ii},2);
  end
end
end
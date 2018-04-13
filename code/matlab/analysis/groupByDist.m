function [inter_cohs, num_in_bins, binned_cohs, std_err] = groupByDist(cohs, dist_edges)
if nargin < 2
  dist_edges = linspace(0, hypot(9,9)*400, 6);
end
if ndims(cohs) == 3
  tmp = zeros([size(cohs),1]); % because reshape(cohs, ..., 1) doesn't work
  tmp(:,:,:,1) = cohs;
  cohs = tmp;
  clear tmp;
end
s = size(cohs);
assert(s(2) == 96 && s(3) == 96);
num_trials = s(end);
dist_mat = funcOverVecs(@dist, 1:96, 1:96);
dist_labs = discretize(dist_mat, dist_edges);
% collapse inner dims for easy indexing
cohs = reshape(cohs, [s(1), 96*96, num_trials]);

nbins = numel(dist_edges)-1;
inter_cohs = zeros(s(1),num_trials,nbins);
num_in_bins = zeros(1,nbins);
binned_cohs = cell(1,nbins);
std_err = zeros(size(inter_cohs));
% just take cohs where i > j to avoid duplicates
mask = logical(tril(ones(96)) - eye(96));
for i=1:nbins
  idxs = and(dist_labs == i, mask);
  num_in_bins(i) = sum(idxs(:));
  binned_cohs{i} = cohs(:,idxs(:),:);
end
for ii=1:nbins
  inter_cohs(:,:,ii) = nanmean(binned_cohs{ii},2); % average all electrode pairs per trial
  % any freq is nan, sum over all electrode pairs
  n_not_nan = squeeze(sum(~any(isnan(binned_cohs{ii}),1),2)); % num_trials x 1
  stds = reshape(nanstd(binned_cohs{ii}, 0, 2), [s(1), num_trials]); % freq x num_trials
  std_err(:,:,ii) = 1.96 * stds  ./ sqrt(n_not_nan'); % implicit expansion
end

  function d = dist(x,y)
    d = norm(getLocFromElectrodeIdx(x) - getLocFromElectrodeIdx(y)) * 400;
  end
end
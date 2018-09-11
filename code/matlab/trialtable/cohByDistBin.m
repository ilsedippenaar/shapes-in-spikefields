function cohs = cohByDistBin(lfps, electrode_mapping, analysis_idxs, dist_edges, mt_params)
% this is for reducing memory load - special case for only getting one coh
% measurement for each distance bin without standard deviation
% lfps = {num_trials} x time x chans
lfps = single(cat(3, lfps{:})); % time x chans x num_trials
lfps = lfps(analysis_idxs, :, :);
if nargin < 4 || isempty(dist_edges)
  dist_idxs = getCohDistIdxs;
else
  dist_idxs = getCohDistIdxs(dist_edges);
end
if nargin < 5 || isempty(mt_params)
  mt_params = getDefaultMtParams;
end
mt_params.trialave = true;
nbins = numel(dist_idxs);

cohs = zeros(numel(getFourierFreqs(analysis_idxs, 1)), nbins);
for i=1:nbins
  lfp_idx_pairs = elecIdxPairToLfpIdxPair(dist_idxs{i}, electrode_mapping);
  for j=1:size(lfp+idx_pairs,1)
    cohs(:, j) = cohs(:,j) + ...
      
  end
end
end
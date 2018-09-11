function cohs = cohByDay(lfps, electrode_mapping, analysis_idxs, mt_params, freq_idxs, idx_pairs)
if nargin < 4 || isempty(mt_params)
  mt_params = getDefaultMtParams;
end
if nargin < 6
  idx_pairs = getCohDistIdxs([0,inf]); % get all non-duplicate electrode pairs
  idx_pairs = idx_pairs{1};
end
% lfps = time x electrode x trial or 1 x trial -> time x electrode
if iscell(lfps)
  lfps = cat(3,lfps{:});
end
if ~isfield(mt_params, 'trialave')
  mt_params.trialave = false;
end

[lfp_idx_pairs, idx_pairs] = elecIdxPairToLfpIdxPair(idx_pairs, electrode_mapping);

num_trials = size(lfps,3);
lfps = single(lfps);

C = coherencyc(lfps(analysis_idxs,1,1), lfps(analysis_idxs,1,1), mt_params);
if nargin < 5 || isempty(freq_idxs)
  freq_idxs = 1:numel(C);
end
if mt_params.trialave
  cohs = nan(numel(freq_idxs), 96, 96);
  for i=1:size(idx_pairs,1)
    elec_i = idx_pairs(i,1);
    elec_j = idx_pairs(i,2);
    coh = coherencyc(lfps(analysis_idxs,lfp_idx_pairs(i,1),:), ...
                     lfps(analysis_idxs,lfp_idx_pairs(i,2),:), mt_params);
    cohs(:,elec_i,elec_j) = coh(freq_idxs);
    cohs(:,elec_j,elec_i) = cohs(:,elec_i,elec_j);
  end
else
  cohs = nan(numel(freq_idxs), num_trials, 96, 96); 
  for i=1:size(idx_pairs,1)
    elec_i = idx_pairs(i,1);
    elec_j = idx_pairs(i,2);
    coh = coherencyc(lfps(analysis_idxs,lfp_idx_pairs(i,1),:), ...
                     lfps(analysis_idxs,lfp_idx_pairs(i,2),:), mt_params);
    cohs(:,:,elec_i,elec_j) = coh(freq_idxs,:);
    cohs(:,:,elec_j,elec_i) = cohs(:,:,elec_i,elec_j);
  end
  cohs = permute(cohs, [1,3,4,2]);
end
end
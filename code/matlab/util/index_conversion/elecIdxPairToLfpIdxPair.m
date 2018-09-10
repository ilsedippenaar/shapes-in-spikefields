function [lfp_idx_pairs, idx_pairs] = elecIdxPairToLfpIdxPair(idx_pairs, electrode_mapping)
% get valid LFP and electrode index pairs
lfp_idxs = arrayfun(@(x) electrode_mapping{x,1}, 1:96, 'UniformOutput', false);
empty_idxs = cellfun(@isempty, lfp_idxs);
lfp_idxs(empty_idxs) = repelem({nan},sum(empty_idxs));
lfp_idxs = [lfp_idxs{:}];
% convert index pairs to lfp index pairs
lfp_idx_pairs = arrayfun(@(x) lfp_idxs(x), idx_pairs);
valid = ~any(isnan(lfp_idx_pairs),2);
idx_pairs = idx_pairs(valid,:);
lfp_idx_pairs = lfp_idx_pairs(valid,:);
end
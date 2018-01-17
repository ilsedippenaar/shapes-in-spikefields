function out = regressReactionTimes(dh, between)
% deprecated
lfps = dh.select('type','lfp', 'between', between, ...
  'trial_section', 'shape_to_saccade', 'trial_result', {'true_positive'});
first_not_empty = find(cellfun(@numel, lfps{1}),1);
last_not_empty = find(cellfun(@numel, lfps{1}), 1, 'last');
out = zeros(last_not_empty-first_not_empty+1, 2, numel(lfps));
for i=1:numel(lfps)
  idx = 1;
  for j=first_not_empty:last_not_empty
    rxn_time = numel(lfps{i}{j});
    t = (1:rxn_time)';
    coeffs = [ones(numel(t),1), t] \ single(lfps{i}{j});
    out(idx,:,i) = [coeffs(2), rxn_time];
    idx = idx + 1;
  end
end
end
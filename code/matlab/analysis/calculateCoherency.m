function cohs = calculateCoherency(lfps, electrode_mappings, mt_params, analysis_window_idxs)
% lfps = days -> electrodes -> time x trials
expanded_lfps = lfps;
for day=1:numel(lfps)
  expanded_lfps{day} = expandCellArray(electrode_mappings{day}, lfps{day});
end
combined_lfps = combineCellArrays('single', expanded_lfps{:});
combined_indices = getCombinedIndices(expanded_lfps);

if nargin < 4
  analysis_window_idxs = 1:size(combined_lfps{1},1);
end

idx = find(~cellfun(@isempty, combined_lfps),1);
C = coherencyc(combined_lfps{idx}(analysis_window_idxs,1), combined_lfps{idx}(analysis_window_idxs,1), mt_params);
cohs = zeros(numel(C), 96, 96);
for i=1:96
  fprintf('%d / %d\n', i, 96);
  cohs(:,i,i) = 1;
  for j=i+1:96
    if isempty(combined_lfps{i}) || isempty(combined_lfps{j})
      C = nan(size(cohs,1),1);
    else
      [~,idxs1,idxs2] = intersect(combined_indices{i}, combined_indices{j});
      C = coherencyc(combined_lfps{i}(analysis_window_idxs,idxs1), ...
                     combined_lfps{j}(analysis_window_idxs,idxs2), mt_params);
    end
    cohs(:,i,j) = C;
    cohs(:,j,i) = C;
  end
end
end
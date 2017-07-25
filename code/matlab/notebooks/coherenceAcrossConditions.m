% for cond=conditions
%   plts = plotCoherenceGraphs(dh, cond{1}, start_time, between(2), 'num_fft', 500, 'window_width', 128, 'freq_bin', [0 25]);
%   fprintf('Saving..\n');
%   saveFigures(plts, fullfile(plot_save_dir, sprintf('coherence_%s.pdf', cond{1})));  
% end

for i=1:numel(conditions)
  cond = conditions{i};
  lfp_slice = dh.getDataSlices(lfps, 'lfp', select_range, cond, 'start_time', start_time);
  n = dh.num_lfp_electrodes;
  plts = gobjects(1,n);
  for j=1:n
    if mod(j,2)
      fprintf('%d...', j);
    end
    if ~mod(j,10)
      fprintf('\n');
    end
    plts(j) = figure('Visible', 'off');
    plotMeanAndStds(lfp_slice{j}, 'x', select_range(1):select_range(2)-1);
    title(sprintf('Electrode %d, N=%d Condition=%s', getElectrodeIdxFromLfpIdx(dh, j), size(lfp_slice{j},2), names{i}));
    plot([-1 -1], ylim(), 'r-');
  end
  saveFigures(plts, fullfile(plot_save_dir, sprintf('traces_1000_%s.pdf', names{i})));
  fprintf('\n');
end
fprintf('All traces saved\n');

for i=1:numel(conditions)
  cond = conditions{i};
  lfp_slices = dh.getDataSlices(lfps, 'lfp', select_range, cond, 'start_time', start_time, 'data_type', 'single');
  [cohs, freqs] = calculateCoherence(lfp_slices, cache_dir, cond, 'num_fft', 500, 'window_width', 128);
  plts = plotCoherenceVariance(cohs, freqs);
  saveFigures(plts, fullfile(plot_save_dir, sprintf('all_std_1000_%s.pdf', names{i})));
end
%% PSD Heatmap with frequency bins
for i=1:numel(lfp_slices)
  fprintf('Plotting PSD heatmaps for condition: %s\n', names{i});
  [psd, f, ~] = calculatePsd(lfp_slices{i}, cache_dir, conditions{i}, 'method', 'mtm');
  plts = gobjects(1,5);
  for j=0:4
    fprintf('%d out of %d\n', j+1, 5);
    plts(j+1) = plotPsdHeatmap(dh, psd, f, [j j+1]*10);
  end
  saveFigures(plts, fullfile(plot_save_dir, sprintf('psd_heatmap_%s.pdf', names{i})));
end

%% Interelectrode distance coherence with distance bins
breaks = linspace(0, hypot(9,9)*400, 6);
dist_bins = [breaks(1:end-1);breaks(2:end)];
for i=1:numel(lfp_slices)
  fprintf('Plotting interelectrode distance coherences for condition: %s\n', names{i});
  [inter_cohs, freqs] = calculateInterelecDistCoherence(dh, lfp_slices{i}, dist_bins, ...
    cache_dir, conditions{i}, 'num_fft', 500, 'window_width', 128);
  plt = plotInterelecDistCoherence(inter_cohs, freqs, dist_bins, 'freq_bin', [0 100]);
  saveFigures(plt, fullfile(plot_save_dir, sprintf('interelec_dist_coh_%s.pdf', names{i})));
end
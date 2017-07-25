%% Initialize
noise_select_params = {'type', 'lfp', 'between', [start_time, size(dh.lfps,1)-30*1000], ...
  'trial_section', 4, 'trial_result', {'true_positive', 'false_negative'}};
shape_select_params = {'type', 'lfp', 'between', [start_time, size(dh.lfps,1)-30*1000], ...
  'trial_section', 5, 'trial_result', {'true_positive', 'false_negative'}};

post_noise_lfps = dh.select(noise_select_params{:});
post_shape_lfps = dh.select(shape_select_params{:});

lfps = post_noise_lfps;
lfp_params = noise_select_params;
name = 'noise';

first_not_empty = find(cellfun(@numel, lfps{1}),1);
last_not_empty = find(cellfun(@numel, lfps{1}),1, 'last');
for i=1:numel(lfps)
  lfps{i} = lfps{i}(first_not_empty:last_not_empty);
end
%% Traces
plts = gobjects(1,numel(lfps));
for i=1:numel(plts)
  plts(i) = plotFullTraces(lfps,i);
end
saveFigures(plts, fullfile(plot_save_dir, 'traces', sprintf('trace_post_%s.pdf', name)));
%% PSD 
[psds, freqs] = calculatePsd(lfps, cache_dir, lfp_params, 'method', 'mtm');
psds = 10*log10(psds);
plt = plotMeanAndStds(psds, 'x', freqs);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
std_plot = figure('Visible', 'off');
plot(freqs, std(psds,0,2) / sqrt(size(psds,2)-1));
title(sprintf('SEM of PSD, N=%d', size(psds,2)));
xlabel('Frequency (Hz)');
ylabel('SEM');
saveFigures([plt, std_plot], fullfile(plot_save_dir, 'psd', sprintf('psd_post_%s.pdf', name)));
%% PSD Heatmap
plts = gobjects(1,10);
for i=1:numel(plts)
  plts(i) = plotPsdHeatmap(dh, psds, freqs, [i-1,i]*5);
end
saveFigures(plts, fullfile(plot_save_dir, 'psd_heatmap', sprintf('psd_heatmap_post_%s.pdf', name)));
%% Coherence and std
[cohs, freqs] = calculateCoherence(lfps, cache_dir, lfp_params, 'method', 'mtm');
plts = plotCoherenceVariance(cohs, freqs);
saveFigures(plts, fullfile(plot_save_dir, 'coherence_var', sprintf('coh_var_post_%s.pdf', name)));
%% Interelectrode distance and coherence
breaks = linspace(0, hypot(9,9)*400, 6);
dist_bins = [breaks(1:end-1);breaks(2:end)];
[inter_cohs, freqs, num_in_bins] = calculateInterelecDistCoherence(dh, lfps, dist_bins, ...
  cache_dir, lfp_params, 'method', 'mtm');
plt = plotInterelecDistCoherence(inter_cohs, freqs, num_in_bins, dist_bins, 'freq_bin', [0 100]);
saveFigures(plt, fullfile(plot_save_dir, 'interelectrode_distance', sprintf('interelec_dist_coh_post_%s.pdf', name)));

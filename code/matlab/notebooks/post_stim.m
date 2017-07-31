%% Initialize
name = 'shape';
if strcmp(name, 'noise')
  lfp_params = {'type', 'lfp', 'between', [start_time, size(dh.lfps,1)-30*1000], ...
  'trial_section', 4, 'trial_result', {'true_positive', 'false_negative'}};
else
  lfp_params = {'type', 'lfp', 'between', [start_time, size(dh.lfps,1)-30*1000], ...
  'trial_section', 5, 'trial_result', {'true_positive'}}; % make sure there is a saccade afterwards
end

lfps = dh.select(lfp_params{:});
first_not_empty = find(cellfun(@numel, lfps{1}),1);
last_not_empty = find(cellfun(@numel, lfps{1}),1, 'last');
for i=1:numel(lfps)
  lfps{i} = lfps{i}(first_not_empty:last_not_empty);
end
%% Traces
plts = gobjects(1,numel(lfps));
for i=1:numel(plts)
  plts(i) = plotFullTraces(lfps{i},i);
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


%% From now on, we deal with the combined LFPs for the given condition (i.e. plotFullTraces)
combined_lfps = cellArray2mat(cellfun(@combineVariableLengthLfps, lfps, 'UniformOutput', false));
%% PSD
mt_params = [];
mt_params.Fs = dh.lfp_sample_freq;
mt_params.pad = 1;
mt_params.err = [1, 0.05]; % theoretical errors at p=0.05
mt_params.tapers = [6,5];
mt_params.trialave = true;
[s,f,Serr] = mtspectrumc(combined_lfps, mt_params);
plt = figure('Visible', 'off');
cutoff_freq = 200;
plot(f(f < cutoff_freq), 10*log10(s(f < cutoff_freq)));
hold on
plot(f(f<cutoff_freq), 10*log10(Serr(:,f<cutoff_freq)'),'r-')
saveFigures(plt, fullfile(plot_save_dir, 'psd', sprintf('post_combined_%s.pdf', name)));
%% Spectrogram
mt_params = [];
mt_params.Fs = dh.lfp_sample_freq;
mt_params.pad = 1;
mt_params.tapers = [6,5];
mt_params.trialave = true;
window_size = 0.064;
[S,t,f] = mtspecgramc(combined_lfps, [window_size, window_size/2], mt_params);
S = S(:,f < cutoff_freq);
plt = plotSpectrogram(S,t,f(f<cutoff_freq));
saveFigures(plt, fullfile(plot_save_dir, 'spectrogram', sprintf('post_%s.pdf', name)));
%% Cosine distance adjacency matrix
% This looks for groupings of electrodes (and possibly how time-domain
% similarities change between the noise and shape stimulus.
dist_mat = calculateAdjMatrix(combined_lfps, 'cosine');
plt = figure('Visible', 'off');
imagesc(dist_mat);
colorbar;
colormap hot;
saveFigures(plt, fullfile(plot_save_dir, 'adjacency_matrix', sprintf('cosine_dist_post_%s.pdf', name)));
%% Factoring the adjacency matrix to look for electrode groupings
[V,d] = eig(dist_mat, 'vector');
n = 20;
per_row = 5;
plt = figure('Visible', 'off');
for i=1:n
  ax = subplot(ceil(n/per_row),per_row,i);
  axis(ax, 'square');
  electrodeHeatmap(electrodeVecToMat(dh,V(:,end-i+1)), ax);
  copyobj(ax,plt);
end
saveFigures(plt, fullfile(plot_save_dir, 'adjacency_matrix_decomposition', ...
  sprintf('cosine_dist_decomp_post_%s.pdf', name)));

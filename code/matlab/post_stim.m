noise_select_params = {'type', 'lfp', 'between', [start_time, size(dh.lfps,1)-30*1000], ...
  'trial_section', 4, 'trial_result', {'true_positive', 'false_negative'}};
shape_select_params = {'type', 'lfp', 'between', [start_time, size(dh.lfps,1)-30*1000], ...
  'trial_section', 5, 'trial_result', {'true_positive', 'false_negative'}};

post_noise_lfps = dh.select(noise_select_params{:});
post_shape_lfps = dh.select(shape_select_params{:});

first_not_empty = find(cellfun(@numel, post_noise_lfps{1}),1);
last_not_empty = find(cellfun(@numel, post_noise_lfps{1}),1, 'last');
for i=1:numel(post_noise_lfps)
  post_noise_lfps{i} = post_noise_lfps{i}(first_not_empty:last_not_empty);
end

first_not_empty = find(cellfun(@numel, post_shape_lfps{1}),1);
last_not_empty = find(cellfun(@numel, post_shape_lfps{1}),1, 'last');
for i=1:numel(post_shape_lfps)
  post_shape_lfps{i} = post_shape_lfps{i}(first_not_empty:last_not_empty);
end
%% Traces
figure(plotFullTraces(post_noise_lfps, 1));
figure(plotFullTraces(post_shape_lfps, 1));
%% Interelectrode distance and coherence
breaks = linspace(0, hypot(9,9)*400, 6);
dist_bins = [breaks(1:end-1);breaks(2:end)];
[inter_cohs, freqs] = calculateInterelecDistCoherence(dh, post_shape_lfps, dist_bins, ...
  cache_dir, shape_select_params, 'method', 'mtm');
plt = plotInterelecDistCoherence(inter_cohs, freqs, dist_bins, 'freq_bin', [0 100]);
figure(plt)
%% Coherence and std

%% PSD heatamps
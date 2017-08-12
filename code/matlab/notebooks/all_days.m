% TODOs:
%   exclude shape data where noise_to_shape is not long enough
%   look at mtspectrumc_unequal_trial_lengths

% Deals only with data from trial sections (i.e. dh.select())
% Make: 
%   traces
%   psds and variance (for each electrode)
%   psd heatmap
%   coherence and variance
%   interelectrode distance coherence
%   spectrograms (especially with centered shape)
%   reaction times versus slope after shape
%   ? time normalized spectrogram
%% Parameters
% Selecting
noise_offset = 1001; % start 1 sec after noise to avoid big onset stimulus
min_trial_lengths = [noise_offset+50, 50]; % at least 50 points for each trial
length_presaccade_response = 150;
% Multitaper
T = 0.5;
W = 6;
% Frequency resolution
num_freqs = 256; % only for (strictly) positive frequencies
% PSD heatmap
freq_bin_size = 10;
% Plotting
freq_cutoff = 120;

%% Initialize
trial_sections = {'noise_to_shape', 'shape_to_saccade'};
select_params = cell(1,2);
select_params{1} = {'type', 'lfp', 'trial_result', {'true_positive'}, 'trial_section', trial_sections{1}};
select_params{2} = select_params{1};
select_params{2}{end} = trial_sections{2};

cache_params = [];
cache_params.data_type = [];
cache_params.dates = arrayfun(@(x) sprintf('%s (%d)', x.date, x.number_on_date), dhs, 'UniformOutput', false);
cache_params.betweens = betweens;
cache_params(2) = cache_params(1);
[cache_params.select_params] = deal(select_params{:});
[cache_params.min_trial_lengths] = deal(min_trial_lengths(1), min_trial_lengths(2));
cache_params(1).noise_offset = noise_offset;

electrode_mapping = cell(96,2);
electrode_mapping(:,1) = num2cell(1:96);

lfps = cell(numel(dhs),numel(trial_sections));
noise_idx = find(strcmp(trial_sections, 'noise_to_shape'));
days_to_delete = false(1,numel(dhs));
for i=1:numel(dhs)
  day_lfps = dhs(i).select('type', 'lfp', 'trial_result', {'true_positive'}, 'trial_section', trial_sections);
  for j=1:numel(day_lfps)
    sections_to_delete = false(1,size(day_lfps{j},1));
    for k=1:size(day_lfps{j},1)
      if numel(day_lfps{j}{k,noise_idx}) < min_trial_lengths(noise_idx)
        sections_to_delete(k) = true;
      else
        day_lfps{j}{k,noise_idx} = day_lfps{j}{k,noise_idx}(noise_offset:end);
      end
    end
    day_lfps{j}(sections_to_delete,:) = [];
  end
  if all(cellfun(@isempty, day_lfps))
    days_to_delete(i) = true;
  else
    for j=1:numel(trial_sections)
      section_lfps = cellfun(@(c) c(:,j), day_lfps, 'UniformOutput', false);
      lfps{i,j} = expandCellArray(dhs(i).electrode_mapping, section_lfps);
    end
  end
end
lfps(days_to_delete,:) = [];
day_idxs = find(~days_to_delete);

% lfps = cell(numel(dhs),numel(trial_sections));
% num_trials_per_electrode = zeros(96,numel(trial_sections));
% max_section_lengths = zeros(size(lfps));
% for i=1:numel(trial_sections)
%   for j=1:numel(dhs)
%     % select data and expand to 96 electrodes
%     lfps{j,i} = expandCellArray(dhs(j).electrode_mapping, dhs(j).select(select_params{i}{:}));
%     % eliminate too small sections
%     trial_sizes = cellfun(@numel, lfps{j,i}{find(~cellfun(@isempty, lfps{j,i}),1)});
%     for k=1:96
%       if ~isempty(lfps{j,i}{k})
%         lfps{j,i}{k}(trial_sizes < min_trial_lengths(i)) = [];
%       end
%     end
%     % eliminate the noise stimulus onset response
%     if strcmp(trial_sections{i}, 'noise_to_shape')
%       lfps{j,i} = funcInCells(lfps{j,i}, @(c) c(noise_offset:end));
%     end
%     % set some variables for convenience
%     num_trials_per_electrode(:,i) = num_trials_per_electrode(:,i) + ...
%       cellfun(@numel, lfps{j,i})';
%     first_not_empty = find(~cellfun(@isempty, lfps{j,i}),1);
%     max_section_lengths(j,i) = max(cellfun(@numel, lfps{j,i}{first_not_empty}));
%   end
% end

% make combined lfps for time analysis
combined_lfps = cell(1,numel(trial_sections));
for i=1:numel(trial_sections)
  combined_lfps{i} = combineCellArrays([], lfps{:,i});
end

% make padded lfps for freqency analysis
num_fft = 2 .^ nextpow2(max(max_section_lengths));
padded_lfps = cell(size(lfps)); % need to keep days separate for coherence
for i=1:numel(trial_sections)
  padded_lfps(:,i) = funcInCells(lfps(:,i), @padData, [], {num_fft(i), 'single'});
  for j=1:size(padded_lfps,1)
    for k=1:96
      padded_lfps{j,i}{k} = cellArray2mat(padded_lfps{j,i}{k}', 'single', 0);
    end
  end
end

% get shpe-centered lfps for spectrogram
shape_centered_lfps = cell(1,numel(dhs));
% make *really* sure nothing else is happening
shape_centered_condition = cell2struct([...
  {'shape', false, [0,1]}; ...
  {'noise', true, [-(256+noise_offset), 256+length_presaccade_response]}; ...
  {'saccade', true, [-(256+noise_offset),256+length_presaccade_response]}], {'vec','negate','range'}, 2);
for i=1:numel(dhs)
  shape_centered_lfps{i} = expandCellArray(dhs(i).electrode_mapping, ...
    dhs(i).getDataSlices('lfp', [-256,256], shape_centered_condition));
end
shape_centered_lfps = combineCellArrays('single', shape_centered_lfps{:});

% get saccade-centered lfps to look at pre-saccadic responses
saccade_centered_lfps = cell(1,numel(dhs));
saccade_centered_condition = cell2struct(...
  {'saccade', false, [0,1]}, {'vec','negate','range'}, 2);
for i=1:numel(dhs)
  saccade_centered_lfps{i} = expandCellArray(dhs(i).electrode_mapping, ...
    dhs(i).getDataSlices('lfp', [-1000,1000], saccade_centered_condition));
end
saccade_centered_lfps = combineCellArrays('single', saccade_centered_lfps{:});
%% Traces
% plot combined traces
for i=1:numel(trial_sections)
  plts = gobjects(1,96);
  for j=1:96
    plts(i) = plotFullTraces(combined_lfps{i}{j}, j);
  end
  saveFigures(plts, fullfile(plot_save_dir, 'all_days/traces', sprintf('traces_%s.pdf', trial_sections{i})));
end
%% PSD
psd_cache_params = cache_params;
[psd_cache_params.date_type] = deal('psd');
all_psds = cell(1,numel(trial_sections));
for i=1:numel(trial_sections)
  combined_padded_lfps = combineCellArrays('single', padded_lfps{:,i});
  [psds, freqs] = calculatePsd(combined_padded_lfps, cache_dir, psd_cache_params(i), ...
    'method', 'mtm', 'T', T, 'W', W);
  factor = (size(psds,1)-1) / num_freqs;
  psds = downsample(psds, factor); % this is necessary so that conditions can be compared
  freqs = downsample(freqs, factor);
  all_psds{i} = psds;
  
  psds = psds(freqs < freq_cutoff,:);
  freqs = freqs(freqs < freq_cutoff);
  psds = 10*log10(psds);
  plt = plotMeanAndStds(psds, 'x', freqs);
  xlabel('Frequency (Hz)');
  ylabel('Power (dB)');
  std_plot = figure('Visible', 'off');
  plot(freqs, std(psds,0,2) / sqrt(size(psds,2)-1));
  title(sprintf('SEM of PSD, N=%d', size(psds,2)));
  xlabel('Frequency (Hz)');
  ylabel('SEM');
  saveFigures([plt, std_plot], fullfile(plot_save_dir, 'all_days/psd', sprintf('psd_%s.pdf', trial_sections{i})));
end
%% PSD heatmap
plts = gobjects(1,5);
range = 10*log10([min(cellfun(@(c) min(min(c)), all_psds)), ...
  max(cellfun(@(c) max(max(c)), all_psds))]);
for i=1:numel(trial_sections)
  for j=1:numel(plts)
    idxs = and(freqs >= (j-1)*freq_bin_size, freqs < j*freq_bin_size);
    psd_map = electrodeVecToMat(electrode_mapping, 10*log10(mean(all_psds{i}(idxs,:),1)));
    plts(j) = electrodeHeatmap(psd_map, [], range);
    title(sprintf('%d to %d Hz', (j-1)*freq_bin_size, j*freq_bin_size));
  end
  saveFigures(plts, fullfile(plot_save_dir, 'all_days/psd_heatmap', sprintf('psd_heatmap_%s.pdf', trial_sections{i})));
end
%% Coherence
coh_cache_params = cache_params;
[coh_cache_params.data_type] = deal('coherence');
all_cohs = cell(1,2);
for i=1:numel(trial_sections)
  num_points = size(padded_lfps{1,i}{1},1)/2+1;
  all_cohs{i} = zeros(num_points,96,96);
  num_counted = zeros(1,96,96); % bookkeeping for weighted average across days
  for j=1:size(padded_lfps,1)
    day_idx = day_idxs(j);
    fprintf('Date: %s (%d)\n', dhs(day_idx).date, dhs(day_idx).number_on_date);
    coh_cache_params(i).date = sprintf('%s (%d)', dhs(day_idx).date, dhs(day_idx).number_on_date);
    [cohs, freqs] = calculateCoherence(padded_lfps{j,i}, cache_dir, coh_cache_params(i), ...
      'method', 'mtm', 'T', T, 'W', W);
    
    first_not_empty = find(~cellfun(@isempty, padded_lfps{j,i}),1);
    num_trials = size(padded_lfps{j,i}{first_not_empty},2);
    all_cohs{i} = all_cohs{i} + cohs*num_trials;
    num_counted = num_counted + (~all(cohs==0,1))*num_trials;
  end
  factor = (size(all_cohs{i},1)-1) / num_freqs;
  all_cohs{i} = downsample( all_cohs{i} ./ num_counted, factor); 
  freqs = downsample(freqs, factor);
  
  idxs = freqs < freq_cutoff;
  plts = plotCoherenceVariance(all_cohs{i}(idxs,:,:),freqs(idxs));
  saveFigures(plts, fullfile(plot_save_dir, 'all_days/coherence', sprintf('coh_%s.pdf', trial_sections{i})));
end
%% Interelectrode coherence
breaks = linspace(0, hypot(9,9)*400, 6);
dist_bins = [breaks(1:end-1);breaks(2:end)];
for i=1:numel(trial_sections)
  [inter_cohs, num_in_bins] = calculateInterelecDistCoherence(all_cohs{i}, dist_bins, electrode_mapping);
  plt = plotInterelecDistCoherence(inter_cohs, freqs, num_in_bins, dist_bins, 'freq_bin', [0 freq_cutoff]);
  saveFigures(plt, fullfile(plot_save_dir, 'all_days/interelectrode_distance', sprintf('interelec_dist_coh_%s.pdf', trial_sections{i})));
end
%% Spectrogram
% plot shape-centered lfps to make sure selection occurred properly
plt = plotMeanAndStds(cellArray2mat(shape_centered_lfps,'int16'));
title(sprintf('N = %d', sum(cellfun(@(c) size(c,2), shape_centered_lfps))));
saveFigures(plt, fullfile(plot_save_dir, 'all_days/traces', 'shape_centered_lfp.pdf'));
% plot spectrogrom
plts = gobjects(1,numel(shape_centered_lfps));
window_size = 0.128;
mt_params = [];
mt_params.tapers = [T*W, 2*T*W-1];
mt_params.fpass = [0,freq_cutoff];
mt_params.Fs = 1000;
mt_params.trialave = true;
for i=1:numel(shape_centered_lfps)
  [S,t,f] = mtspecgramc(shape_centered_lfps{i},[window_size, window_size/2],mt_params);
  plts(i) = plotSpectrogram(S,t,f);
  title(sprintf('Spectrogram for electrode %d', i));
end
saveFigures(plts, fullfile(plot_save_dir, 'all_days/spectrogram', 'shape_centered_spectrogram.pdf'));
%% Reaction times
% plot saccade-cenetered lfps to see when pre-saccadic response occurs
plt = plotMeanAndStds(cellArray2mat(saccade_centered_lfps,'int16'));
title(sprintf('N = %d', sum(cellfun(@(c) size(c,2), saccade_centered_lfps))));
saveFigures(plt, fullfile(plot_save_dir, 'all_days/traces', 'saccade_centered_lfp.pdf'));
% make csv file with cols electrode number, slope, and reaction time to
% plotting in R
shape_idx = find(strcmp(trial_sections, 'shape_to_saccade'));
rxns_and_slopes = cell(1,96);
for i=1:96
  valid_idxs = find(cellfun(@numel, combined_lfps{shape_idx}{i}) >= length_presaccade_response+100); % at least 100 ms before saccade response to regress with
  rxns_and_slopes{i} = zeros(numel(valid_idxs),2);
  for j=1:numel(valid_idxs)
    idx = valid_idxs(j);
    rxn_time = numel(combined_lfps{shape_idx}{i}{idx});
    t = (1:numel(combined_lfps{shape_idx}{i}{idx})-length_presaccade_response)';
    coeffs = [ones(numel(t),1), t] \ single(combined_lfps{shape_idx}{i}{idx}(1:end-length_presaccade_response));
    rxns_and_slopes{i}(j,:) = [coeffs(2), rxn_time];
  end
end
csvwrite(fullfile(data_dir, '../R/', 'rxns_and_slopes.csv'), ...
  [repelem((1:96)', cellfun(@(c) size(c,1), rxns_and_slopes)), ...
  cellArray2mat(rxns_and_slopes')]);
%% Time-normalized spectrogram
plts = gobjects(1,96);
shape_idx = find(strcmp(trial_sections, 'shape_to_saccade'));
window_size = 0.128;
mt_params = [];
mt_params.tapers = [T*W, 2*T*W-1];
mt_params.fpass = [0,freq_cutoff];
mt_params.Fs = 1000;
mt_params.trialave = false;
for i=1:96
  all_specgrams = cell(1,numel(combined_lfps{shape_idx}{i}));
  for j=1:numel(combined_lfps{shape_idx}{i})
    [S,~,f] = mtspecgramc(combined_lfps{i}{j},[window_size, window_size/2],mt_params);
    all_specgrams{j} = S;
  end
  normalized_specgram = normalizeData(all_specgrams);
  t = linspace(0,1,size(normalized_specgram,1));
  plts(i) = plotSpectrogram(S,t,f);
  title(sprintf('Normalized spectrogram for electrode %d', i));
end
saveFigures(plts, fullfile(plot_save_dir, 'all_days/spectrogram', 'shape_to_saccade_normalized_spectrogram.pdf'));
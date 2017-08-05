% Deals only with data from trial sections (i.e. dh.select())
% Make: 
%   psds and variance (for each electrode)
%   psd heatmap
%   interelectrode distance coherence and variance
%   spectrograms (especially with centered shape)
%   reaction times versus slope after shape
%% Parameters
% Selecting
min_trial_lengths = [450, 50];
noise_offset = 401; % start 400 ms after noise to avoid big onset stimulus
% Multitaper
T = 0.5;
W = 6;

%% Initialize
trial_sections = {'noise_to_shape', 'shape_to_saccade'};
select_params = cell(1,2);
select_params{1} = {'type', 'lfp', 'trial_result', {'true_positive'}, 'trial_section', trial_sections{1}};
select_params{2} = select_params{1};
select_params{2}{end} = trial_sections{2};

cache_params = [];
cache_params.dates = arrayfun(@(x) sprintf('%s (%d)', x.date, x.number_on_date), dhs, 'UniformOutput', false);
cache_params.betweens = betweens;
cache_params(2) = cache_params(1);
cache_params(1).select_params = select_params{1};
cache_params(2).select_params = select_params{2};

lfps = cell(numel(dhs),numel(trial_sections));
num_trials_per_electrode = zeros(96,numel(trial_sections));
max_section_lengths = zeros(size(lfps));
for i=1:numel(trial_sections)
  for j=1:numel(dhs)
    % select data and expand to 96 electrodes
    lfps{j,i} = expandCellArray(dhs(j), dhs(j).select(select_params{i}{:}));
    % eliminate too small sections
    trial_sizes = cellfun(@numel, lfps{j,i}{find(~cellfun(@isempty, lfps{j,i}),1)});
    for k=1:96
      if ~isempty(lfps{j,i}{k})
        lfps{j,i}{k}(trial_sizes < min_trial_lengths(i)) = [];
      end
    end
    % eliminate the first 400 ms after noise stimulus to avoid large onset
    % response
    if strcmp(trial_sections{i}, 'noise_to_shape')
      lfps{j,i} = funcInCells(lfps{j,i}, @(c) c(noise_offset:end));
    end
    % set some variables for convenience
    num_trials_per_electrode(:,i) = num_trials_per_electrode(:,i) + ...
      cellfun(@numel, lfps{j,i})';
    first_not_empty = find(~cellfun(@isempty, lfps{j,i}),1);
    max_section_lengths(j,i) = max(cellfun(@numel, lfps{j,i}{first_not_empty}));
  end
end

% make padded lfps for freqency analysis
num_fft = 2 .^ nextpow2(max(max_section_lengths));
padded_lfps = cell(size(lfps));
for i=1:numel(trial_sections)
  padded_lfps(:,i) = funcInCells(lfps(:,i), @padData, [], {num_fft(i), 'single'});
  for j=1:numel(dhs)
    for k=1:96
      padded_lfps{j,i}{k} = cellArray2mat(padded_lfps{j,i}{k}', 'single', 0);
    end
  end
end

clear trial_sizes first_not_empty
%% PSD
all_psds = cell(1,numel(trial_sections));
for i=1:numel(trial_sections)
  combined_lfps = combineCellArrays('single', padded_lfps{:,i});
  [psds, freqs] = calculatePsd(combined_lfps, cache_dir, cache_params, ...
    'method', 'mtm', 'T', T, 'W', W);
  all_psds{i} = psds;
  
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
clear curr_section psds mt_params combined_lfps
%% PSD heatmap
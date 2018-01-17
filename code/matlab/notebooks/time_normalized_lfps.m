% For true positive, false positive, and false negative trials results, plot 
%  - shape- and saccade-centered LFPs
%  - time-normalized trace for shape-to-saccade LFP
%  - time-normalized spectrogram for shape-to-saccade LFP
%  - time-normalized (interelectrode-distance) coheregram for shape-to-saccade LFP
% Look at differences before and after shape stimulus given that the
%   pre-shape signal is stationary and the pre-saccade response has not
%   occurred
% Analyze "covariance corrected" shape-to-saccade signal

%% Parameters
% Selecting
length_postnoise_response = 1001; % start 1 sec after noise to avoid big onset stimulus
length_presaccade_response = 150;
% Time normalization
normalize_length = 1e3;
% Frequency analysis and multitaper
T = 0.5;
W = 6;
freq_cutoff = 120;
% Mock electrode mapping for identity case (i.e. expanded cell arrays)
electrode_mapping = cell(96,2);
electrode_mapping(:,1) = num2cell(1:96);
%% Initialize
fprintf('Intializing...\n');
trial_results = {'true_positive', 'false_positive', 'false_negative'};
align_names = {'shape', 'saccade'};

trial_result_vecs = cell(numel(dhs), numel(align_names),numel(trial_results));
all_valid_trials = cell(size(trial_result_vecs));
for i=1:numel(trial_results)
  for j=1:numel(align_names)
    for k=1:numel(dhs)
      valid_trials = strcmp({dhs(k).trials.result}, trial_results{i});
      if strcmp(trial_results{i}, 'false_positive')
        % some shapes occur after the saccade, so make sure those trials
        % are excluded
        all_events = [dhs(k).trials.sections];
        valid_trials = and(valid_trials, ...
          cellfun(@isempty, all_events(dhs(k).section_map('shape_to_stop'):dhs(k).num_trial_sections+1:end)));
      end
      if strcmp(align_names{j}, 'saccade')
        align_events = {dhs(k).trials.saccade};
      else
        all_events = [dhs(k).trials.sections];
        align_events = all_events(dhs(k).section_map('shape_to_stop'):dhs(k).num_trial_sections+1:end);
      end
      trial_result_vecs{k,j,i} = [align_events{valid_trials}];
      all_valid_trials{k,j,i} = find(valid_trials);
    end
  end
end

names = {[trial_results{1}, '_', align_names{1}], ...
         [trial_results{1}, '_', align_names{2}], ...
         [trial_results{2}, '_', align_names{2}], ...
         [trial_results{3}, '_', align_names{1}]};
trial_result_vecs = [trial_result_vecs(:,1,1), ...
                     trial_result_vecs(:,2,1), ...
                     trial_result_vecs(:,2,2), ...
                     trial_result_vecs(:,1,3)];
all_valid_trials = [all_valid_trials(:,1,1), ...
                     all_valid_trials(:,2,1), ...
                     all_valid_trials(:,2,2), ...
                     all_valid_trials(:,1,3)];

% Get shape-to-saccade lfps for each day and combine them
fprintf('Selecting shape to saccade lfps...\n');
day_lfps = cell(1,size(all_valid_trials,1));
for i=1:size(all_valid_trials,1)
  fprintf('%d / %d\n', i, size(all_valid_trials,1));
  day_lfps{i} = expandCellArray(dhs(i).electrode_mapping, ...
    dhs(i).select('type', 'lfp', 'trial_nums', all_valid_trials{i,1}, 'trial_section', 'shape_to_saccade'));
end
s_to_s_lfps = combineCellArrays([], day_lfps{:});

combined_indices = cell(1,96);
num_in_day = cellfun(@(c) max(cellfun(@numel, c)), day_lfps);
for i=1:96
  combined_indices{i} = zeros(1,sum(num_in_day));
  curr_idx = 1;
  for j=1:numel(day_lfps)
    if ~isempty(day_lfps{j}{i})
      combined_indices{i}(curr_idx:curr_idx+num_in_day(j)-1) = 1;
    end
    curr_idx = curr_idx + num_in_day(j);
  end
  combined_indices{i} = find(combined_indices{i});
end

% Normalize lfps
fprintf('Normalizing LFPs...\n');
normalized_lfps = cell(1,96);
for i=1:96
  normalized_lfps{i} = cellArray2mat(normalizeData(s_to_s_lfps{i},'len', normalize_length, 'method', 'linear')')';
end
%% Aligned LFPs
fprintf('Plotting aligned LFPs for each condition...\n');
for i=1:size(trial_result_vecs,2)
  fprintf('Aligning and plotting %s LFPs\n', names{i});
  lfps = cell(1,size(trial_result_vecs,1));
  conditions = cell2struct([ ...
    {[], false, [0,1]}; ...
    {'noise', true, [-length_postnoise_response-250,0]}], {'vec','negate','range'}, 2);
  for j=1:size(trial_result_vecs,1)
    fprintf('%d / %d\n', j, size(trial_result_vecs,1));
    conditions(1).vec = trial_result_vecs{j,i};
    lfps{j} = expandCellArray(dhs(j).electrode_mapping, ...
      dhs(j).getDataSlices('lfp', [-250,250], conditions));
  end
  lfps = combineCellArrays('int16', lfps{:});
  num_counted = max(cellfun(@(c) size(c,2), lfps));
  lfps = cellArray2mat(lfps, 'int16');
  plt = plotMeanAndStds(lfps, 'x', -250:249);
  title(sprintf('N = %d', num_counted));
  xlabel('Time (ms)');
  saveFigures(plt, fullfile(plot_save_dir, 'all_days/traces/aligned', sprintf('%s_centered_lfp.pdf', names{i})));
end
%% Time-normalized traces
fprintf('Plotting time-normalized traces...\n');
plts = gobjects(1,96);
for i=1:96
  fprintf('%d / %d\n', i, 96);
  plts(i) = plotMeanAndStds(normalized_lfps{i}, 'x', linspace(0,1,normalize_length));
  title(sprintf('Electrode %d', i));
end
saveFigures(plts, fullfile(plot_save_dir, 'all_days/traces', 'shape_to_saccade_normalized_trace.pdf'));
%% Time-normalized spectrogram
fprintf('Plotting time-normalized spectrogram...\n');
mt_params = [];
mt_params.tapers = [T*W, 2*T*W-1];
mt_params.fpass = [0,freq_cutoff];
mt_params.Fs = normalize_length;
mt_params.trialave = true;
window_size = 1/20; % 1/20th of the total size
plts = gobjects(1,96);
for i=1:96
  fprintf('%d / %d\n', i, 96);
  [S,t,f] = mtspecgramc(normalized_lfps{i}, [window_size, window_size/2],mt_params);
  plts(i) = plotSpectrogram(S,t,f);
  title(sprintf('Electrode %d', i));
end
saveFigures(plts, fullfile(plot_save_dir, 'all_days/spectrogram', 'shape_to_saccade_normalized_spectrogram2.pdf'));
%% Time-normalized coheregram
fprintf('Plotting time-normalized coheregram...\n');
window_size = 1/20; % 1/20th of the total size
mt_params = [];
mt_params.tapers = dpsschk([T*W, 2*T*W-1], normalize_length*window_size, normalize_length);
mt_params.fpass = [0,freq_cutoff];
mt_params.Fs = normalize_length;
mt_params.trialave = true;

% Calculate coheregrams
[~,~,~,~,~,t,f] = cohgramc(mean(normalized_lfps{1},2), mean(normalized_lfps{2},2),  [window_size, window_size/2], mt_params);
all_cohgrams = zeros(numel(t),numel(f),96,96);
for i=1:96
  fprintf('%d / %d\n', i, 96);
  all_cohgrams(:,:,i,i) = 1;
  for j=i+1:96
    [~,idxs1,idxs2] = intersect(combined_indices{i},combined_indices{j});
    [C,~,~,~,~,t,f] = cohgramc(normalized_lfps{i}(:,idxs1), normalized_lfps{j}(:,idxs2),  [window_size, window_size/2], mt_params);
    all_cohgrams(:,:,j,i) = C;
    all_cohgrams(:,:,i,j) = C;
  end
end

% Separate into distance bins
breaks = linspace(0, hypot(9,9)*400, 6);
dist_bins = [breaks(1:end-1);breaks(2:end)];
[inter_cohgrams, num_in_bins] = calculateInterelecDistCoherence(all_cohgrams, dist_bins, electrode_mapping);

% Plot and save
plts = gobjects(1,size(inter_cohgrams,3));
for i=1:size(inter_cohgrams,3)
  plts(i) = plotSpectrogram(inter_cohgrams(:,:,i),t,f,false);
  title(sprintf('%.0f - %.0f $\\mu m$, n=%d', dist_bins(1,i), dist_bins(2,i), num_in_bins(i)), 'Interpreter', 'latex');
end
saveFigures(plts, fullfile(plot_save_dir, 'all_days/cohgram', 'shape_to_saccade_time_normalized_cohgram.pdf'));
%% Interelectrode distance coherence with normalized data
window_size = 0.064;
mt_params = [];
mt_params.tapers = [T*W, 2*T*W-1];
mt_params.fpass = [0,freq_cutoff];
mt_params.Fs = 1000;
mt_params.trialave = false;
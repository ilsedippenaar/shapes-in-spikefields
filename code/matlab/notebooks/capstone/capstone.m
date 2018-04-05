%% Parameters
folder_name = 'capstone';
freq_cutoff = 100;
alpha_range = [7.5,12.5];
beta_range = [12.6,30];
gamma_range = [31,100];
all_ranges = [alpha_range; beta_range; gamma_range];

monkey_name = dhs(1).monkey_name;

names = {'success', 'failure'};
select_range = [-256, 256];
conds = {
  cell2struct([...
    {'shape', false, [0,1]}; ...
    {'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]}; ...
    {'saccade', true, [-1000, params.length_min_reaction_time]}; ... % shapes can come on *after* a saccade too
    {'saccade', false, [params.length_min_reaction_time, 1000]}], {'vec','negate','range'}, 2), ...
  cell2struct([...
    {'shape', false, [0,1]}; ...
    {'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]}; ...
    {'saccade', true, [-1000, 1000]}], {'vec','negate','range'}, 2)};

mt_params = [];
mt_params.tapers = [params.T*params.W, 2*params.T*params.W-1];
mt_params.fpass = [0,params.freq_cutoff];
mt_params.Fs = params.Fs;
mt_params.trialave = true;
%% Initialize
grouped_dhs = groupDhs(dhs);
electrode_mappings = cellfun(@(c) c(1).electrode_mapping, grouped_dhs, 'UniformOutput', false);
electrode_mappings = arrayfun(@(c) electrode_mappings, 1:numel(conds), 'UniformOutput', false);

fprintf('Selecting data...\n');
spikes = cell(size(conds));
lfps = cell(size(conds));
for cond_num=1:numel(conds)
  spikes{cond_num} = cell(1,numel(grouped_dhs));
  lfps{cond_num} = cell(1,numel(grouped_dhs));
  for day=1:numel(grouped_dhs)
    num_in_day = numel(grouped_dhs{day});
    fprintf('%s, day %d / %d\n', names{cond_num}, day, numel(grouped_dhs));
    day_lfps = cell(1,num_in_day);
    day_spikes = cell(1,num_in_day);
    for num=1:num_in_day
      day_lfps{num} =   grouped_dhs{day}(num).getDataSlices('lfp', select_range, conds{cond_num});
      day_spikes{num} = grouped_dhs{day}(num).getDataSlices('spike', select_range, conds{cond_num});
    end
    lfps{cond_num}{day} =   combineCellArrays('single', day_lfps{:});
    spikes{cond_num}{day} = combineCellArrays([], day_spikes{:});
  end
  valid_days = ~cellfun(@(c) isempty(c) || all(cellfun(@isempty, c)), spikes{cond_num});
  spikes{cond_num} = spikes{cond_num}(valid_days);
  lfps{cond_num} = lfps{cond_num}(valid_days);
  electrode_mappings{cond_num} = electrode_mappings{cond_num}(valid_days);
end

fprintf('\nDetermining trials with relatively high LFP power...\n');
analysis_window_idxs = 1:-select_range(1);
% calculate average PSD for each trial
psds = cell(size(lfps));
for i=1:numel(lfps) % result
  psds{i} = cell(size(lfps{i}));
  for j=1:numel(lfps{i}) % day
    fprintf('%d\n', j);
    % rearrange to have electrodes nested in a trial
    trial_lfps = cell(1,numel(lfps{i}{j}));
    for k=1:numel(lfps{i}{j}) % electrode
      % trials to cell array
      trial_lfps{k} = mat2cell(lfps{i}{j}{k}, size(lfps{i}{j}{k},1), repelem(1,size(lfps{i}{j}{k},2)));
    end
    trial_lfps = combineCellArrays('single', trial_lfps{:});
    psds{i}{j} = cell(1,numel(trial_lfps));
    for k=1:numel(trial_lfps)
      % look from the beginning of the data up to when the shape comes on
      [psds{i}{j}{k},pre_shape_freqs] = mtspectrumc(trial_lfps{k}(analysis_window_idxs,:), mt_params);
      psds{i}{j}{k} = psds{i}{j}{k} / norm(psds{i}{j}{k}, Inf); % normalize so max = 1
    end
    psds{i}{j} = double([psds{i}{j}{:}]);
  end
end

% label trials by the relative prevalence of a given frequency
num_splits = 3;
combined_psds = cell(size(psds));
for i=1:numel(combined_psds)
  combined_psds{i} = [psds{i}{:}];
end
trial_pow_labs = psds; % take advantage of copy-on-write semantics
lfp_pow_bins = quantile([combined_psds{:}],(0:num_splits)/num_splits, 2);
for i=1:numel(conds)
  for j=1:numel(psds{i}) % days
    for k=1:size(psds{i}{j},1) % freqs
      trial_pow_labs{i}{j}(k,:) = discretize(psds{i}{j}(k,:), lfp_pow_bins(k,:));
    end
  end
end

% label trials by the relative prevalence in a frequency range
range_psds = cell(size(combined_psds));
for i=1:numel(range_psds)
  range_psds{i} = cell(size(all_ranges,1),1);
  for j=1:size(all_ranges,1)
    idxs = pre_shape_freqs >= all_ranges(j,1) & pre_shape_freqs < all_ranges(j,2);
    range_psds{i}{j} = mean(combined_psds{1}(idxs,:),1);
  end
  range_psds{i} = vertcat(range_psds{i}{:});
end
trial_pow_labs_range = cell(1,numel(conds));
lfp_pow_bins = quantile([range_psds{:}],(0:num_splits)/num_splits, 2);
for i=1:numel(conds)
  trial_pow_labs_range{i} = cell(1,numel(psds{i}));
  for j=1:numel(psds{i}) % days
    trial_pow_labs_range{i}{j} = zeros(size(all_ranges,1),size(psds{i}{j},2));
    for k=1:size(all_ranges,1) % freq ranges
      idxs = pre_shape_freqs >= all_ranges(k,1) & pre_shape_freqs < all_ranges(k,2);
      trial_pow_labs_range{i}{j}(k,:) = discretize(mean(psds{i}{j}(idxs,:),1), lfp_pow_bins(k,:));
    end
  end
end

fprintf('\nCalculating coherences...\n');
cohs = cell(size(conds));
for cond_num=1:numel(conds)
  cohs{cond_num} = calculateCoherency(lfps{cond_num}, electrode_mappings{cond_num}, mt_params, analysis_window_idxs);
end

fprintf('Calculating phases at shape onset...\n');
all_freq_ranges = [alpha_range; beta_range; gamma_range]';
phases = cell(2,3); % sucess / failure by alpha/beta/gamma
for i=1:3
  [~, split_lfps] = splitByPowLabels(spikes, lfps, trial_pow_labs_range, num_splits, i); % 1 = alpha, 2 = beta, 3 = gamma
  mean_pre_shape_freqs = conv(pre_shape_freqs, [1,1]/2, 'valid');
  idxs = and(all_freq_ranges(1,i) < mean_pre_shape_freqs, mean_pre_shape_freqs < all_freq_ranges(2,i));
  for j=1:2
    % {3,1} = high power, success, {3,2} = high power, failure
    phases{j,i} = getLfpPhases(split_lfps{3,j}, pre_shape_freqs, -select_range(1)+1, [1,-select_range(1)], params.Fs, 'fft');
    phases{j,i} = angle([phases{j,i}{idxs}]); % combine fft results in bin of interest
  end
end

fprintf('\nCalculating spatial distribution of coherences...\n');
coh_sigmas = cellfun(@getSpatialCohErrors, cohs, 'UniformOutput', false);

fprintf('Calculating STAs...\n');
timwin = [-0.1, 0.1];
stas = cell(2,2); % success / failure by alpha / gamma
pow_idxs = [1,3;
            1,3];
success_idxs = [1,1;
                2,2];
for i=1:numel(stas)
  [split_spikes, split_lfps] = splitByPowLabels(spikes, lfps, trial_pow_labs_range, num_splits, pow_idxs(i));
  [pow_spikes, pow_lfps] = deal(split_spikes{3,success_idxs(i)}, split_lfps{3,success_idxs(i)}); % 3 = high power
  [y, x, std_err] = calculateSTA(pow_lfps, pow_spikes, select_range(1), params.Fs, timwin);
  stas{i} = {x,y,std_err};
end

fprintf('Calculating PPC values for all units in all conditions...\n');
ppcs = cell(2,2); % same as STA
pow_idxs = [1,3;
            1,3];
success_idxs = [1,1;
                2,2];
for i=1:numel(ppcs)
  [split_spikes, split_lfps] = splitByPowLabels(spikes, lfps, trial_pow_labs_range, num_splits, pow_idxs(i));
  [pow_spikes, pow_lfps] = deal(split_spikes{3,success_idxs(i)}, split_lfps{3,success_idxs(i)}); % 3 = high power
  [ppc, f] = ...
    calculatePPC(pow_lfps, pow_spikes, electrode_mappings{success_idxs(i)}, select_range(1), params.Fs);
  ppcs{i} = {f, ppc};
end

fprintf('Making trial table...\n');
trial_tbl = arrayfun(@(x) x.toTrialTable, dhs, 'UniformOutput', false)';
trial_tbl = vertcat(trial_tbl{:});
trial_tbl = [trial_tbl table((1:size(trial_tbl,1))', 'VariableNames', {'abs_num'})];

% there is a shape
% the trial has enough data to select (probably unnecessary)
% shape comes on late enough
% there isn't a saccade too early
valid_trials = [~isnan(trial_tbl.shape), ...
                all(trial_tbl.start < trial_tbl.shape + select_range < trial_tbl.stop,2), ... 
                ~between(trial_tbl.noise, trial_tbl.shape-params.length_postnoise_response+select_range(1), trial_tbl.shape+select_range(2)), ...
                ~between(trial_tbl.saccade , trial_tbl.shape-1000, trial_tbl.shape+params.length_min_reaction_time)];
valid_trials = all(valid_trials, 2);
analysis_tbl = trial_tbl(valid_trials,:);

% select data
idxs = select_range(1):select_range(2)-1;
analysis_lfps = rowfun(@(x, start, shape) {x{1}(idxs+shape-start+1,:)}, ...
  analysis_tbl, 'InputVariables', {'lfps', 'start', 'shape'});
analysis_tbl.lfps = table2cell(analysis_lfps);

select_spikes = @(x, shape, range) x(binarySearch(x, range(1), ']', true):binarySearch(x, range(2), '(', true))-shape;
analysis_spikes = rowfun(@(x, shape) funcInCells(x, select_spikes, [], {shape,idxs([1,end])+shape}), ...
  analysis_tbl, 'InputVariables', {'spikes', 'shape'});
analysis_tbl.spikes = table2cell(analysis_spikes);
%% Plot 1 - data summary
if strcmp(monkey_name, 'jaws')
  % A - trial structure summary (success)
  [events, trial_lfp, trial_spikes] = getSummaryData(dhs(3), 2);
  saveToR([], [], [], data_dir, folder_name, '1_a_trial', 'trial', events);
  saveToR(1:numel(trial_lfp), trial_lfp, [], data_dir, folder_name, '1_a_lfp');
  saveToR(1:numel(trial_lfp), trial_spikes, [], data_dir, folder_name, '1_a_spikes');

  % B - trial structure summary (failure)
  [events, trial_lfp, trial_spikes] = getSummaryData(dhs(3), 19);
  saveToR([], [], [], data_dir, folder_name, '1_b_trial', 'trial', events);
  saveToR(1:numel(trial_lfp), trial_lfp, [], data_dir, folder_name, '1_b_lfp');
  saveToR(1:numel(trial_lfp), trial_spikes, [], data_dir, folder_name, '1_b_spikes');
end

% C - firing rate distribution
rates = cell(1,numel(dhs));
for i=1:numel(dhs)
  num_fired = cellfun(@numel, dhs(i).spikes);
  rates{i} = num_fired ./ size(dhs(i).lfps,1) * 1000; % num spikes / second
end
rates = [rates{:}];
saveToR([],[],[], data_dir, folder_name, sprintf('1_c_%s', monkey_name), 'rates', rates);

% D - reaction time distribution
rxn_times = cell(1,numel(dhs));
for i=1:numel(dhs)
  rxn_times{i} = cell(1,dhs(i).num_trials);
  for j=1:dhs(i).num_trials
    rxn_times{i}{j} = dhs(i).trials(j).saccade - dhs(i).trials(j).sections{5};
  end
  rxn_times{i} = [rxn_times{i}{:}];
end
rxn_times = [rxn_times{:}];
saveToR([], [], [], data_dir, folder_name, sprintf('1_d_%s', monkey_name), 'rxn_time', rxn_times);
%% Plot 2 - info analysis summary

%% Plot 3 - LFPs and frequency
% A - performance and presence of frequency
n_success = size([trial_pow_labs{1}{:}],2);
n_failure = size([trial_pow_labs{2}{:}],2);
success_counts = countcats(categorical([trial_pow_labs{1}{:}]),2);
failure_counts = countcats(categorical([trial_pow_labs{2}{:}]),2);
norm_success_counts = success_counts / n_success;
norm_failure_counts = failure_counts / n_failure;
success_rate = norm_success_counts ./ (norm_success_counts+norm_failure_counts);
x = pre_shape_freqs(pre_shape_freqs < freq_cutoff);
y = success_rate(pre_shape_freqs < freq_cutoff,:);
% take binomial (bernoulli) error for success rates
n = success_counts + failure_counts;
std_err = sqrt(success_rate .* (1-success_rate) ./ n);
saveToR(x, y, std_err, data_dir, folder_name, sprintf('3_a_%s', monkey_name));

% B - phase histograms
saveToR(phases, [], [], data_dir, folder_name, sprintf('3_b_%s', monkey_name));

% C - interelectrode distance coherences
breaks = linspace(0, hypot(9,9)*400, 6);
dist_bins = [breaks(1:end-1);breaks(2:end)];
electrode_mapping = cell(96,2);
electrode_mapping(:,1) = num2cell(1:96);
% coherences are averaged across trials, so we can take a weighted average
% over trials to get the coherence considered across all conditions
n_success = sum(cellfun(@(d) size(d{1},2), lfps{1}));
n_failure = sum(cellfun(@(d) size(d{1},2), lfps{2}));
all_cohs = sum(cat(4, cohs{1}*n_success, cohs{2}*n_failure),4) / (n_success+n_failure);
[y, ~, ~, std_err] = calculateInterelecDistCoherence(all_cohs, dist_bins, electrode_mapping);
x = pre_shape_freqs(pre_shape_freqs < freq_cutoff);
saveToR(x, y, std_err, data_dir, folder_name, sprintf('3_c_%s', monkey_name), 'dist_bins', dist_bins);

% D - sptial coherence and performance
x = pre_shape_freqs(pre_shape_freqs < freq_cutoff);
y = [coh_sigmas{1}(:,1), coh_sigmas{2}(:,2)];
std_err = [coh_sigmas{1}(:,2:3), coh_sigmas{2}(:,2:3)];
saveToR(x, y, std_err, data_dir, folder_name, sprintf('3_d_%s', monkey_name));
%% Plot 4 - spike triggered averages
letters = {'a','b','c','d'};
for i=1:numel(letters)
  [x,y,std_err] = deal(stas{i}{:});
  saveToR(x, y, std_err, data_dir, folder_name, sprintf('4_%s_%s', letters{i}, monkey_name));
end
%% Plot 5 - PPCs
letters = {'a','b','c','d'};
for i=1:numel(letters)
  [x,y] = deal(ppcs{i}{:});
  x = mean([x{:}],2);
  all_ppcs = [y{:}];
  all_ppcs = all_ppcs(:,~any(isnan(all_ppcs),1));
  y = mean(all_ppcs,2);
  std_err = 1.96 * std(all_ppcs, 0, 2) / sqrt(size(all_ppcs,2));
  saveToR(x, y, std_err, data_dir, folder_name, sprintf('5_%s_%s', letters{i}, monkey_name));
end
%% Plot 6 - learning
% A - success rate over time
ma_win = 100; % moving average
success = contains(analysis_tbl.result, 'true');
x = analysis_tbl.date(floor(movmean((1:size(analysis_tbl,1))', ma_win)));
y = movmean(success, ma_win);
saveToR(x, y, [], data_dir, folder_name, sprintf('6_a_%s', monkey_name));

% B - spatial coherence over time
sigma_by_row = @(x) {getSpatialCohErrors(calculateCoherencyFromRows(x(:,1),x(:,2),mt_params,analysis_window_idxs))};
row_sigmas = splitapply(sigma_by_row, ...
  [analysis_tbl.lfps, analysis_tbl.electrode_mapping], ...
  findgroups(analysis_tbl.date));

x = unique(analysis_tbl.date);
freqs = pre_shape_freqs(pre_shape_freqs < freq_cutoff);
y = cat(3, row_sigmas{:});
y = squeeze(y(:,1,:)); % just get the mean sigmas, not the std errs
saveToR(x, y, [], data_dir, folder_name, sprintf('6_b_%s', monkey_name), 'freqs', freqs);

% C - PPC over time
[day_ppcs, ppc_freqs] = splitapply(@(x) calculatePPCFromRows(x(:,1), x(:,2), x(:,3), select_range(1), params.Fs), ...
  [analysis_tbl.lfps, analysis_tbl.spikes, analysis_tbl.electrode_mapping], ...
  findgroups(analysis_tbl.date));

x = unique(analysis_tbl.date);
freqs = ppc_freqs{1}(:,1);
y = funcInCells(day_ppcs, @(x) nanmean(x,2), @() nan(numel(freqs),1));
y = [y{:}];
saveToR(x, y, [], data_dir, folder_name, sprintf('6_c_%s', monkey_name), 'freqs', freqs);

% D - variance / std of PPC over time
x = unique(analysis_tbl.date);
freqs = ppc_freqs{1}(:,1);
y = funcInCells(day_ppcs, @(x) nanstd(x,0,2), @() nan(numel(freqs),1));
y = [y{:}];
saveToR(x, y, [], data_dir, folder_name, sprintf('6_d_%s', monkey_name), 'freqs', freqs);
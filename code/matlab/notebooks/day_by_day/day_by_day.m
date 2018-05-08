monkey_name = 'jaws';
analysis_type = 'shape';
folder_name = 'day_by_day';

alpha_range = [7.5,25];
gamma_range = [31,100];
all_ranges = [alpha_range; gamma_range];

freq_cutoff = 100;

% sorting into power bins
num_splits = 3;

% moving average params
ma_window_size = 300;
ma_step_size = 50; 

% coherence distance params
dist_edges = linspace(0, hypot(9,9)*400, 6);

trial_dir = fullfile(data_dir, 'trials', monkey_name);
lfp_dir = fullfile(data_dir, 'lfps', monkey_name);
dh_dir = fullfile(cache_dir, 'data_handlers', monkey_name);
save_dir = fullfile(cache_dir, 'trial_tables', analysis_type, monkey_name);
all_save_dir = fullfile(cache_dir, 'trial_tables', 'all', monkey_name);
output_dir = fullfile(cache_dir, 'output', monkey_name);

mt_params = [];
mt_params.tapers = [params.T*params.W, 2*params.T*params.W-1];
mt_params.Fs = params.Fs;
mt_params.trialave = true;

if strcmp(monkey_name, 'jaws')
  dates = datetime('2013-05-1'):datetime('2013-07-19');
  learn_start = datetime('2013-06-3');
  %dates = datetime('2013-06-3'):datetime('2013-07-19');
else
  dates = datetime('2011-7-22'):datetime('2011-12-20');
  learn_start = datetime('2011-11-7');
  %dates = datetime('2011-11-7'):datetime('2011-12-20');
end

if strcmp(analysis_type, 'shape')
  analysis_idxs = 1:256;
else
  analysis_idxs = 1:256;
end
%%
tbl = makeTrialTable(dates, trial_dir, lfp_dir, dh_dir, save_dir, analysis_type, ...
  'clean', true, ...
  'min_reaction_time', params.length_min_reaction_time, ...
  'monkey_name', monkey_name);
grouped = groupBy(tbl, 'date');
%% Add analysis columns for all data
% Success
tbl.success = strcmp(tbl.result, 'true_positive');

% Reaction time
tbl.rxn_time = tbl.saccade - tbl.shape;

% PSD
tbl.psd = getPsd(tbl, analysis_idxs, mt_params);

% Power bins (all)
tbl.pow_labs = getLabels(tbl, 'psd', num_splits);

% Power bin (low and high)
freqs = getFourierFreqs(analysis_idxs, params.Fs);
binned_pows = zeros(size(tbl,1), size(all_ranges,1));
for i=1:size(all_ranges,1)
  idxs = between(freqs, all_ranges(i,1), all_ranges(i,2));
  binned_pows(:,i) = mean(tbl.psd(:,idxs),2);
end
pow_bins = quantile(binned_pows, (0:num_splits)/num_splits, 1);
pow_labs = apply(@(i) discretize(binned_pows(:,i), pow_bins(:,i)), ...
                 1:size(all_ranges,1), 1, size(tbl,1)); % discretize each row (row = trial)
tbl.low_high_pow_labs = pow_labs;

% LFP phases
% taking the mean of a bunch of phases - usually > 70 so, the resultant
% length (mod of the circular mean) isn't expected to be 0 (mean approx.
% 0.106) - 95% quantile is approx. 0.207
phases = rowfun(@(x) {fft(single(x{1}(analysis_idxs,:)))}, tbl, 'InputVariables', {'lfps'});
phases = rowfun(@(x) {x{1}(1:floor(end/2)+1,:)}, phases);
tbl.phases = phases.Var1;

% Coherence (mean and stderr in far bin)
% grouping by day necessary for doing coherence with one consistent electrode mapping
% don't average trials
trial_cohs = cell(1,numel(grouped));
trial_mt_params = mt_params;
trial_mt_params.trialave = false;
% get farthest distance bin
%idx_pairs = getCohDistIdxs(dist_edges);
%idx_pairs = idx_pairs{end};
parfor ii=1:numel(grouped)
  fprintf("Coherence: day %s\n", grouped{ii}(1,:).date{1});
  % coh = freq x 96 x 96 x num_trials
  coh = cohByDay(grouped{ii}.lfps, grouped{ii}(1,:).electrode_mapping{1}, ...
    analysis_idxs, trial_mt_params);
  binned_cohs = groupCohByDist(coh, dist_edges); % {dist_bin} x freqs x pairs x trials
  trial_cohs{ii} = cell(1,numel(binned_cohs));
  for j=1:numel(binned_cohs)
    [m,s] = calcMeanAndStd(binned_cohs{j},2);
    trial_cohs{ii}{j} = [m,s]; % freqs x 2 x trials (horzcat works on dim 2, which we averaged over)
  end
  trial_cohs{ii} = cat(4,trial_cohs{ii}{:}); % freqs x 2 x trials x dist_bin
  trial_cohs{ii} = permute(trial_cohs{ii}, [1,2,4,3]);
end
trial_cohs = cat(4,trial_cohs{:}); % freqs x 2 x dist_bin x num trials
tbl.coh = squeeze(permute(trial_cohs(:,1,:,:), [4,1,2,3]));
tbl.coh_std_err = squeeze(permute(trial_cohs(:,2,:,:), [4,1,2,3]));

coh_bins = quantile(tbl.coh(:,:,end), (0:num_splits)/num_splits, 1);
tbl.far_coh_labs = apply(@(i) discretize(tbl.coh(:,i,end), coh_bins(:,i)), ...
                 1:size(tbl.coh,2), 1, size(tbl,1)); % discretize each row (row = trial)

% Informative neurons
% PPC of informative neurons
ppcs = cell(size(tbl,1),1);
parfor ii=1:size(tbl,1)
  ppcs(ii) = calculatePPCFromRows(tbl(ii,:).lfps, tbl(ii,:).spikes, ...
    tbl(ii,:).electrode_mapping, -256, params.Fs);
end
tbl.ppcs = ppcs; % {trials} x freqs x num_units

freqs = 0:10:100;
mean_pl = cellfun(@(c) calcMeanAndStd(c, 2, true), tbl.ppcs, 'UniformOutput', false);
for i=1:numel(mean_pl)
  if size(mean_pl{i},1) ~= numel(freqs)
    mean_pl{i} = nan(numel(freqs),1);
  end
end
mean_pl = [mean_pl{:}]';

tbl.pl_labs = getLabels(cell2struct({mean_pl}, 'pl'), 'pl', num_splits);
%% Plot 2 - Data Summary
if strcmp(monkey_name, 'jaws')
  tmp_trial_dir = fullfile(data_dir, 'trials', monkey_name);
  tmp_lfp_dir = fullfile(data_dir, 'lfps', monkey_name);
  tmp_dh_dir = fullfile(cache_dir, 'data_handlers', monkey_name);
  dhs = DataHandler.fromDates('2013-05-1', tmp_trial_dir, tmp_lfp_dir, tmp_dh_dir, ...
    'clean', true, ...
    'min_reaction_time', params.length_min_reaction_time, ...
    'monkey_name', monkey_name);
  dhs = dhs.split();
  dhs = [dhs{:}];
  dhs = dhs([dhs.num_trials] ~= 0);
  
  % A - trial structure summary (success)
  [events, trial_lfp, trial_spikes] = getSummaryData(dhs(3), 2);
  saveToR([], [], [], data_dir, folder_name, '2_a_trial', 'trial', events);
  saveToR(1:numel(trial_lfp), trial_lfp, [], data_dir, folder_name, '2_a_lfp');
  saveToR(1:numel(trial_lfp), trial_spikes, [], data_dir, folder_name, '2_a_spikes');

  % B - trial structure summary (failure)
  [events, trial_lfp, trial_spikes] = getSummaryData(dhs(3), 19);
  saveToR([], [], [], data_dir, folder_name, '2_b_trial', 'trial', events);
  saveToR(1:numel(trial_lfp), trial_lfp, [], data_dir, folder_name, '2_b_lfp');
  saveToR(1:numel(trial_lfp), trial_spikes, [], data_dir, folder_name, '2_b_spikes');
  clear dhs tmp_trial_dir tmp_lfp_dir tmp_dh_dir
end

% C - reaction times
saveToR([], [], [], data_dir, folder_name, sprintf('2_c_%s', monkey_name), 'rxn_time', tbl.rxn_time);

% D - success rate vs shape - noise time
% need to get all the trials, not just "shape" trials to calculate proper
% stimulus time
all_tbl = makeTrialTable(dates, trial_dir, lfp_dir, dh_dir, all_save_dir, [], ...
  'clean', true, ...
  'min_reaction_time', params.length_min_reaction_time, ...
  'monkey_name', monkey_name);
all_tbl = all_tbl(all_tbl.shape - all_tbl.noise > 0,:);
all_tbl = all_tbl(contains(all_tbl.result, {'true_positive', 'false_negative'}),:);
all_tbl.success = strcmp(all_tbl.result, 'true_positive');

num_bins = 50; % for quantiles
bin_size = 100; % for fixed bin sizes

stim_times = accumarray(all_tbl.success + 1, all_tbl.shape - all_tbl.noise, [], @(x) {x});
stim_times{1} = stim_times{1}(~isnan(stim_times{1})); % if any false negative trial creep in

%edges = quantile(vertcat(stim_times{:}), num_bins);
edges = 0:bin_size:ceil(max(vertcat(stim_times{:}))/bin_size)*bin_size;
counts = cellfun(@(c) histcounts(c, edges), stim_times, 'UniformOutput', false);
counts = vertcat(counts{:});

n = sum(counts,1)';
x = conv(edges, [1,1]/2, 'valid')';
y = counts(2,:)' ./ n;
std_err = sqrt(y .* (1-y) ./ n);

saveToR(x, y, std_err, data_dir, folder_name, sprintf('2_d_%s', monkey_name));
clear all_tbl
%% Plot 3 - Frequencies
% A - Success rate by frequency split by LFP power
[x,y,std_err] = getFreqSuccess(tbl, 'pow_labs', analysis_idxs, freq_cutoff, params.Fs);
saveToR(x, y, std_err, data_dir, folder_name, sprintf('3_a_%s', monkey_name));

% TODO: bootstrap all of the confidence intervals
% idxs = getFourierFreqs(analysis_idxs, params.Fs) < freq_cutoff;
% ci = bootci(round(size(tbl,1)), ...
%   {@getFreqSuccessBootstrap, tbl.success, tbl.pow_labs(:,idxs)}, ...
%   'Options', statset('UseParallel', true), 'type', 'per');
% plot(x,y)
% hold on
% plot(x, squeeze(ci(:,3,:)),'r--')
% 
% x = tbl.psd(tbl.success,26);
% y = tbl.psd(~tbl.success,26);
% x = x(randi(numel(x), 1, numel(y)));
% signrank(x, y)

% calculate chi-squared test statistic
p = zeros(1,11);
for i=1:11
  [~,~,p(i)] = crosstab(tbl.success, tbl.pl_labs(:,i));
end

% B - Phases in bins of interest (low / high)
% TODO: calculate PLV or other phase-locking measure here
grouped_phases = cell(2, size(all_ranges,1));
for i=1:size(all_ranges,1)
  idxs = between(freqs, all_ranges(i,1), all_ranges(i,2));
  high_pow_tbl = tbl(tbl.low_high_pow_labs(:,i)==num_splits,:);
  phases = rowfun(@(x) {angle(mean(x{1}(idxs,:),2))}, high_pow_tbl, 'InputVariables', {'phases'});
  grouped_phases{1,i} = vertcat(phases(high_pow_tbl.success,:).Var1{:});
  grouped_phases{2,i} = vertcat(phases(~high_pow_tbl.success,:).Var1{:});
end
saveToR(grouped_phases, [], [], data_dir, folder_name, sprintf('3_b_%s', monkey_name));

% C - "Confidence" (i.e. reaction time) in power bin
s_tbl = tbl(tbl.success,:);
% sort into groups based on pow_lab
rxn_groups_low = accumarray(s_tbl.low_high_pow_labs(:,1), s_tbl.rxn_time, [], @(x) {x});
rxn_groups_high = accumarray(s_tbl.low_high_pow_labs(:,2), s_tbl.rxn_time, [], @(x) {x});
saveToR(1:num_splits, {rxn_groups_low, rxn_groups_high}, [], ...
        data_dir, folder_name, sprintf('3_c_%s', monkey_name));
%% Plot 4 - Coherence
% A - Coherence binned by distance
x = getFourierFreqs(analysis_idxs, params.Fs);
[y,std_err] = calcMeanAndStd(tbl.coh,1);
y = squeeze(y);
std_err = squeeze(std_err);
idxs = x < freq_cutoff;
x = x(idxs);
y = y(idxs,:);
std_err = std_err(idxs,:);
saveToR(x, y, std_err, data_dir, folder_name, sprintf('4_a_%s', monkey_name), 'dist_edges', dist_edges);

% B - Coherence Distribtuion
% TODO: change coherence distribution calculation


% C - Success rate by frequency of far coherence
[x,y,std_err] = getFreqSuccess(tbl, 'far_coh_labs', analysis_idxs, freq_cutoff, params.Fs);
saveToR(x, y, std_err, data_dir, folder_name, sprintf('4_c_%s', monkey_name));
%% Plot 5: Spike density
win_size = 50;
smoothed = cellfun(@(c) mean(calculateSmoothedSpikes(c, gausswin(win_size), -256, 512),2), tbl.spikes, 'UniformOutput', false);
smoothed = [smoothed{:}];

% A - Spike density in success
[y,std_err] = calcMeanAndStd(smoothed(:,tbl.success),2);
[~,x] = calculateSmoothedSpikes({},gausswin(win_size),-256,512);
saveToR(x, y, std_err, data_dir, folder_name, sprintf('5_a_%s', monkey_name));

% B - Spke density in failure
[y,std_err] = calcMeanAndStd(smoothed(:,~tbl.success),2);
[~,x] = calculateSmoothedSpikes({},gausswin(win_size),-256,512);
saveToR(x, y, std_err, data_dir, folder_name, sprintf('5_b_%s', monkey_name));
%% Plot 6: Phase locking
% A - Phase locking in success / failure
x = getFourierFreqs(1:100, params.Fs);
x = x(x <= freq_cutoff)';
mean_ppcs = cellfun(@(c) calcMeanAndStd(c,2), tbl.ppcs, 'UniformOutput', false);
s_ppcs = mean_ppcs(tbl.success);
f_ppcs = mean_ppcs(~tbl.success);

s_ppcs = s_ppcs(cellfun(@numel, s_ppcs) == numel(x)); % some are empty
s_ppcs = [s_ppcs{:}];
[s_y,s_std_err] = calcMeanAndStd(s_ppcs,2,true);

f_ppcs = f_ppcs(cellfun(@numel, f_ppcs) == numel(x));
f_ppcs = [f_ppcs{:}];
[f_y,f_std_err] = calcMeanAndStd(f_ppcs,2,true);

saveToR(x, [s_y,f_y], [s_std_err,f_std_err], data_dir, folder_name, sprintf('6_a_%s', monkey_name));

apply(@(i) ranksum(s_ppcs(i,:),f_ppcs(i,:)), 1:11, 1, 1)

% B - Success vs phase locking bin
[x,y,std_err] = getFreqSuccess(tbl, 'pl_labs', 1:100, freq_cutoff, params.Fs);
saveToR(x, y, std_err, data_dir, folder_name, sprintf('6_b_%s', monkey_name));
%% Plot 7: Learning
learn_tbl = tbl(and(tbl.shapecoh < 80, ...
                    datetime(tbl.date) >= learn_start), :);

x = movmean(1:size(learn_tbl,1), ma_window_size, 'EndPoints', 'discard');
x = x(1:ma_step_size:end);

group_bounds = [round(x)-ma_window_size/2; round(x)+ma_window_size/2-1];
group_idxs = apply(@(x) x(1):x(2), group_bounds, 1, ma_window_size);
trial_groups = cell(1,size(group_idxs,2));
for i=1:numel(trial_groups)
  trial_groups{i} = learn_tbl(group_idxs(:,i),:);
end

success = learn_tbl.success;
suc = movmean(success, ma_window_size, 'EndPoints', 'discard');
suc = suc(1:ma_step_size:end);

% A - learning rate
%success = all_tbl(all_tbl.shapecoh < 70,:).success;
y = movmean(success, ma_window_size, 'EndPoints', 'discard');
y = y(1:ma_step_size:end);
saveToR(x', y, [], data_dir, folder_name, sprintf('7_a_%s', monkey_name));

% B - frequency change over time (low power portion of plot 3a over time)
trial_group_success = cell(1,numel(trial_groups));
for i=1:numel(trial_group_success)
  [~,trial_group_success{i},~] = getFreqSuccess(trial_groups{i}, 'pow_labs', analysis_idxs, freq_cutoff, params.Fs);
end
trial_group_success = cat(3, trial_group_success{:});

y = squeeze(trial_group_success(:,1,:)); % get the lowest power bin
y = y - repmat(suc, 1, size(y,1))'; % calculate delta
freqs = getFourierFreqs(analysis_idxs, params.Fs);
freqs = freqs(freqs < freq_cutoff);

saveToR(x, y, [], data_dir, folder_name, sprintf('7_b_%s', monkey_name), 'freqs', freqs);

% C - frequency change over time (high power portion of 3a over time)
y = squeeze(trial_group_success(:,end,:)); % get the highest power bin
y = y - repmat(suc, 1, size(y,1))'; % calculate delta
freqs = getFourierFreqs(analysis_idxs, params.Fs);
freqs = freqs(freqs < freq_cutoff);
saveToR(x, y, [], data_dir, folder_name, sprintf('7_c_%s', monkey_name), 'freqs', freqs);
%% Day analyses - require multiple trials

% % Coherence distribution - takes too long with individual trials
% coh_distribution = cell(1,numel(grouped));
% parfor ii=1:numel(grouped)
%   fprintf("Day %s\n", grouped{ii}(1,:).date{1});
%   day_cohs = cohByDay(grouped{ii}.lfps, grouped{ii}(1,:).electrode_mapping{1}, analysis_idxs, mt_params);
%   coh_distribution{ii} = getSpatialCohErrors(day_cohs);
% end
% x = unique(tbl.date);
% freqs = getFourierFreqs(analysis_idxs,1000);
% y = cat(3, coh_distribution{:});
% y = squeeze(y(:,1,:)); % just get the mean sigmas, not the std errs
% saveToR(x, y, [], data_dir, folder_name, sprintf('coh_distrib_day_%s_%s', analysis_type, monkey_name), 'freqs', freqs);
% 
% 
% % Phase locking
% ppcs = cell(1,numel(grouped));
% freqs = cell(size(ppcs));
% parfor ii=1:numel(grouped)
%   [ppcs(ii), freqs(ii)] = calculatePPCFromRows(grouped{ii}.lfps, grouped{ii}.spikes, ...
%     grouped{ii}(1,:).electrode_mapping, -256, params.Fs);
% end
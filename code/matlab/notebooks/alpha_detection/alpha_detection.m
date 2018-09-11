monkey_name = 'jaws';
analysis_type = 'preshape';
folder_name = 'alpha_detection';

alpha_range = [15, 20];
gamma_range = [31,100];
all_ranges = [alpha_range; gamma_range];

freq_cutoff = 100;

% sorting into power bins
num_splits = 2;

% moving average params
ma_window_size = 300;
ma_step_size = 50; 

% coherence distance params
dist_edges = linspace(0, hypot(9,9)*400, 6);

% spike density calculation
win_size = 50;

if strcmp(analysis_type, 'shape')
  analysis_idxs = 1:256;
  trial_start = -256;
  analysis_save_dir = 'shape';
elseif strcmp(analysis_type, 'postshape')
  analysis_idxs = 256:512;
  trial_start = -256;
  analysis_save_dir = 'shape';
elseif strcmp(analysis_type, 'prenoise')
  analysis_idxs = 1:500;
  trial_start = -500;
  analysis_save_dir = 'prenoise';
end

trial_dir = fullfile(data_dir, 'trials', monkey_name);
lfp_dir = fullfile(data_dir, 'lfps', monkey_name);
dh_dir = fullfile(cache_dir, 'data_handlers', monkey_name);
save_dir = fullfile(cache_dir, 'trial_tables', analysis_save_dir, monkey_name);
all_save_dir = fullfile(cache_dir, 'trial_tables', 'all', monkey_name);
analyzed_dir = fullfile(cache_dir, 'trial_tables', 'analyzed', monkey_name);

mt_params = [];
mt_params.tapers = dpss(numel(analysis_idxs), params.T*params.W, 2*params.T*params.W-1);
mt_params.Fs = params.Fs;
mt_params.trialave = true;
mt_params.pad = -1; % don't pad

if strcmp(monkey_name, 'jaws')
  dates = datetime('2013-05-1'):datetime('2013-07-19');
  learn_start = datetime('2013-06-3');
  %dates = datetime('2013-06-3'):datetime('2013-07-19');
else
  dates = datetime('2011-7-22'):datetime('2011-12-20');
  learn_start = datetime('2011-11-7');
  %dates = datetime('2011-11-7'):datetime('2011-12-20');
end
%%
analyzed_name = fullfile(analyzed_dir, sprintf('%s.mat', analysis_type));
if exist(analyzed_name, 'file') ~= 2
  tbl = makeTrialTable(dates, trial_dir, lfp_dir, dh_dir, save_dir, analysis_save_dir, ...
  'clean', true, ...
  'min_reaction_time', params.length_min_reaction_time, ...
  'monkey_name', monkey_name);
  tbl.lfps = cellfun(@(c) mean(c,2), tbl.lfps, 'UniformOutput', false);
  
  dates = unique(datetime(tbl.date));
  for i=1:numel(dates)
    info_units = getVisuallyResponsiveUnits(dates(i), monkey_name, params.length_min_reaction_time, ...
      trial_dir, lfp_dir, dh_dir, all_save_dir);
    idxs = datetime(tbl.date) == dates(i);
    if numel(tbl(idxs,:).spikes{1}) == numel(info_units)
      fprintf('Eliminating %d units for day %s...\n', sum(~info_units), char(dates(i)));
      tbl(idxs,:).spikes = cellfun(@(c) c(info_units), tbl(idxs,:).spikes, 'UniformOutput', false);
    end
  end
  
  grouped = groupBy(tbl, 'date');
  
  % Add analysis columns
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
  
  % PSTHs (technically spike densities)
  total_length = size(tbl.lfps{1},1);
  smoothed = cellfun(@(c) mean(...
    calculateSmoothedSpikes(c, gausswin(win_size), trial_start, total_length), ...
    2) / numel(c), tbl.spikes, 'UniformOutput', false);
  tbl.psth = [smoothed{:}]';

  n = size(tbl,1);
  dist_idxs  = getCohDistIdxs([0, inf]);
  dist_idxs = dist_idxs{1}; % percalculating this saves approx. 1 sec per loop
  trial_cohs = cell(n,1);
  parfor ii=1:n
    fprintf("Coherence: trial %d / %d\n", ii, n);
    % coh = freq x 96 x 96
    coh = cohByDay(tbl(ii,:).lfps, tbl(ii,:).electrode_mapping{1}, ...
      analysis_idxs, mt_params, [], dist_idxs);
    binned_cohs = groupCohByDist(coh, dist_edges); % {dist_bin} x freqs x pairs
    binned_cohs = cellfun(@(c) calcMeanAndStd(c,2), binned_cohs, 'UniformOutput', false);
    trial_cohs{ii} = [binned_cohs{:}]; % freq x dist_bin
  end
  tbl.coh = permute(cat(3,trial_cohs{:}), [3,1,2]); % n x freq x dist_bin

  coh_bins = quantile(tbl.coh(:,:,end), (0:num_splits)/num_splits, 1);
  tbl.far_coh_labs = apply(@(i) discretize(tbl.coh(:,i,end), coh_bins(:,i)), ...
                   1:size(tbl.coh,2), 1, size(tbl,1)); % discretize each row (row = trial)
  
  if strcmp(analysis_type, 'shape')
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
  end
  
  save(analyzed_name, 'tbl', '-v7.3');
else
  tbl = load(analyzed_name);
  tbl = tbl.tbl;
end
%% Plot 1: Coherence
x = getFourierFreqs(analysis_idxs, params.Fs);
[y,std_err] = calcMeanAndStd(tbl.coh,1);
y = squeeze(y);
std_err = squeeze(std_err);
idxs = x < freq_cutoff;
x = x(idxs);
y = y(idxs,:);
std_err = std_err(idxs,:);
saveToR(x, y, std_err, data_dir, folder_name, sprintf('1_%s_%s', analysis_type, monkey_name), 'dist_edges', dist_edges);
%% Plot 2: Performance vs LFP alpha power and phase
% A - Success vs LFP power
[x,y,std_err] = getFreqSuccess(tbl, 'pow_labs', analysis_idxs, freq_cutoff, params.Fs);
saveToR(x, y, std_err, data_dir, folder_name, sprintf('2_a_%s_%s', analysis_type, monkey_name));

% B - Alpha phase
% TODO: instead of angle of resultant length of amplitude normalized LFP
% spectrum, calculate some phase locking measure between LFPs
freqs = getFourierFreqs(analysis_idxs, params.Fs);
grouped_phases = cell(2, 1);
idxs = between(freqs, alpha_range(1), alpha_range(2));
% take circular mean of LFP phases
phases = cellfun(@(c) angle(mean(c(idxs,:),2)), tbl.phases, 'UniformOutput', false);
grouped_phases{1} = vertcat(phases{tbl.success}); % combine different frequencies together
grouped_phases{2} = vertcat(phases{~tbl.success});
saveToR(grouped_phases, [], [], data_dir, folder_name, sprintf('2_b_%s_%s', analysis_type, monkey_name));

% C - how different are the phase distributions at a frequency?
freqs = getFourierFreqs(analysis_idxs, params.Fs);
p_vals = zeros(1,numel(freqs));
phases = cellfun(@(c) angle(mean(c, 2)), tbl.phases, 'UniformOutput', false);
phases = [phases{:}];
s_phases = phases(:,tbl.success);
f_phases = phases(:,~tbl.success);
for i=1:numel(p_vals)
  [~,p_vals(i)] = circDistribTest(s_phases(i,:), f_phases(i,:));
end
idxs = freqs < freq_cutoff;
plt = figure;
plot(freqs(idxs), p_vals(idxs));
hold on
plot([0,freqs(find(idxs,1,'last'))], [0.05,0.05], 'b--')
title('Prob that success and failure phase distributions are equal')
xlabel('Frequency (Hz)')
ylabel('p value')
saveFigures(plt, fullfile(plot_save_dir, folder_name, sprintf('2_c_%s_%s.pdf', analysis_type, monkey_name)), false);
%% Plot 3: Spike densities
if strcmp(analysis_type, 'shape')
  total_length = size(tbl.lfps{1},1);
  valid_trials = cellfun(@(c) ~all(cellfun(@(d) isempty(d) || all(d > 0), c)), tbl.spikes);
  valid_trials = and(valid_trials, tbl.low_high_pow_labs(:,1) == num_splits);
  % TODO: look to see how robust PSTH phase relationship is
  valid_trials = and(valid_trials, randi(10,size(tbl,1),1) > 2); % 80%
  spike_tbl = tbl(valid_trials,:);
  
  % A - spike densities in success / failure
  %idxs = true(size(tbl,1),1);
  [s_y,s_std_err] = calcMeanAndStd(spike_tbl.psth(spike_tbl.success, :),1);
  [f_y,f_std_err] = calcMeanAndStd(spike_tbl.psth(~spike_tbl.success,:),1);
  [~,x] = calculateSmoothedSpikes({}, gausswin(win_size), trial_start, total_length);
  saveToR(x', [s_y', f_y'], [s_std_err', f_std_err'], ...
    data_dir, folder_name, sprintf('3_a_%s_%s', analysis_type, monkey_name));

  % B - phases of PSTH
  [~,x] = calculateSmoothedSpikes({}, gausswin(win_size), trial_start, total_length);
  phases = fft(spike_tbl.psth(:,x<0)');
  phases = phases(1:floor(end/2)+1,:);
  
  s = spike_tbl.success;
  high = spike_tbl.low_high_pow_labs(:,1) == num_splits;
  %figure,histogram(angle(phases(10,and(s,high))))
  %figure,histogram(angle(phases(10,and(~s,high))))
  
  freqs = getFourierFreqs(x(x<0), params.Fs);
  p_vals = zeros(1,numel(freqs));
  s_phases = angle(phases(:,and(s,high)));
  f_phases = angle(phases(:,and(~s,high)));
  for i=1:numel(p_vals)
    [~,p_vals(i)] = circDistribTest(s_phases(i,:), f_phases(i,:));
  end
  idxs = freqs < freq_cutoff;
  plt = figure('Visible', 'off');
  plot(freqs(idxs), p_vals(idxs));
  hold on
  plot([0,freqs(find(idxs,1,'last'))], [0.05,0.05], 'b--')
  title('Prob that success and failure phase distributions are equal')
  xlabel('Frequency (Hz)')
  ylabel('p value')
  saveFigures(plt, fullfile(plot_save_dir, folder_name, sprintf('3_b_%s_%s.pdf', analysis_type, monkey_name)));
  
  [~,idx] = min(abs(freqs - 17));
  plt = figure;
  plot(phases(idx,s),'r.')
  hold on
  plot(phases(idx,~s),'b.')
  legend({'Success', 'Failure'})
  xlabel('Real'); ylabel('Imaginary');
  saveFigures(plt, fullfile(plot_save_dir, folder_name, sprintf('3_c_%s_%s.png', analysis_type, monkey_name)));
end
%% Plot 4: Phase locking in low  / high alpha power bin
if any(strcmp(tbl.Properties.VariableNames, 'ppcs'))
  x = getFourierFreqs(1:100, params.Fs);
  x = x(x <= freq_cutoff)';
  mean_ppcs = cellfun(@(c) calcMeanAndStd(c,2), tbl.ppcs, 'UniformOutput', false);

  % PPCs of low alpha power trials
  s_low_pow_ppcs = mean_ppcs(and(tbl.low_high_pow_labs(:,1) == 1, tbl.success));
  f_low_pow_ppcs = mean_ppcs(and(tbl.low_high_pow_labs(:,1) == 1, ~tbl.success));

  s_low_pow_ppcs = s_low_pow_ppcs(cellfun(@numel, s_low_pow_ppcs) == numel(x)); % some are empty
  s_low_pow_ppcs = [s_low_pow_ppcs{:}];
  [s_y,s_std_err] = calcMeanAndStd(s_low_pow_ppcs,2,true);

  f_low_pow_ppcs = f_low_pow_ppcs(cellfun(@numel, f_low_pow_ppcs) == numel(x)); % some are empty
  f_low_pow_ppcs = [f_low_pow_ppcs{:}];
  [f_y,f_std_err] = calcMeanAndStd(f_low_pow_ppcs,2,true);

  saveToR(x,[s_y, f_y], [s_std_err, f_std_err], data_dir, folder_name, sprintf('4_a_%s_%s', analysis_type, monkey_name));

  % PPCs of high alpha trials
  s_high_pow_ppcs = mean_ppcs(and(tbl.low_high_pow_labs(:,1) == 2, tbl.success));
  f_high_pow_ppcs = mean_ppcs(and(tbl.low_high_pow_labs(:,1) == 2, ~tbl.success));

  s_high_pow_ppcs = s_high_pow_ppcs(cellfun(@numel, s_high_pow_ppcs) == numel(x)); % some are empty
  s_high_pow_ppcs = [s_high_pow_ppcs{:}];
  [s_y,s_std_err] = calcMeanAndStd(s_high_pow_ppcs,2,true);

  f_high_pow_ppcs = f_high_pow_ppcs(cellfun(@numel, f_high_pow_ppcs) == numel(x)); % some are empty
  f_high_pow_ppcs = [f_high_pow_ppcs{:}];
  [f_y,f_std_err] = calcMeanAndStd(f_high_pow_ppcs,2,true);

  saveToR(x,[s_y, f_y], [s_std_err, f_std_err], data_dir, folder_name, sprintf('4_b_%s_%s', analysis_type, monkey_name));
end
%% Plot 5: Informative neurons and alpha phase / power
%% Does prenoise gamma power index "eagerness?" i.e. false alarm rate
% can compare saccade to no saccade condition or correct and false alarm
% rate
if strcmp(analysis_type, 'prenoise')
  % false alarm tbl
  fa_tbl = tbl(or(strcmp(tbl.result, 'true_positive'), strcmp(tbl.result, 'false_negative')),:);
  fa_tbl = fa_tbl(datetime(fa_tbl.date) >= learn_start,:);
  [x,y,std_err] = getFreqSuccess(fa_tbl, 'pow_labs', analysis_idxs, freq_cutoff, params.Fs);
  saveToR(x, y, std_err, data_dir, folder_name, sprintf('false_neg_%s_%s', analysis_type, monkey_name));
end
%% Does prenoise gamma power correlate with reaction time?
if strcmp(analysis_type, 'prenoise')
  s_tbl = tbl(tbl.success,:);
  rxn_groups_low = accumarray(s_tbl.low_high_pow_labs(:,1), s_tbl.rxn_time, [], @(x) {x});
  rxn_groups_high = accumarray(s_tbl.low_high_pow_labs(:,2), s_tbl.rxn_time, [], @(x) {x});
  saveToR(1:num_splits, {rxn_groups_low, rxn_groups_high}, [], ...
          data_dir, folder_name, sprintf('reaction_time_%s', monkey_name));
  
  fa_tbl = tbl(or(strcmp(tbl.result, 'true_positive'), strcmp(tbl.result, 'false_positive')),:);
  y = fa_tbl.saccade - fa_tbl.noise;
  rxn_groups_low = accumarray(fa_tbl.low_high_pow_labs(:,1), y, [], @(x) {x});
  rxn_groups_high = accumarray(fa_tbl.low_high_pow_labs(:,2), y, [], @(x) {x});
  saveToR(1:num_splits, {rxn_groups_low, rxn_groups_high}, [], ...
          data_dir, folder_name, sprintf('reaction_time2_%s', monkey_name));
end
%% Does the role of alpha/beta/gamma power change over the course of learning?
if strcmp(analysis_type, 'prenoise')
  fa_tbl = tbl(or(strcmp(tbl.result, 'true_positive'), strcmp(tbl.result, 'false_positive')),:);
  fa_tbl = fa_tbl(and(fa_tbl.shapecoh < 80, datetime(fa_tbl.date) >= learn_start),:);
  
  x = movmean(1:size(fa_tbl,1), ma_window_size, 'EndPoints', 'discard');
  x = x(1:ma_step_size:end);

  group_bounds = [round(x)-ma_window_size/2; round(x)+ma_window_size/2-1];
  group_idxs = apply(@(x) x(1):x(2), group_bounds, 1, ma_window_size);
  trial_groups = cell(1,size(group_idxs,2));
  for i=1:numel(trial_groups)
    trial_groups{i} = fa_tbl(group_idxs(:,i),:);
  end

  success = fa_tbl.success;
  suc = movmean(success, ma_window_size, 'EndPoints', 'discard');
  suc = suc(1:ma_step_size:end);

  % A - learning rate
  y = movmean(success, ma_window_size, 'EndPoints', 'discard');
  y = y(1:ma_step_size:end);
  figure;
  plot(x', y)
  title('Learning rate')
  xlabel('Trial Index'); ylabel('Success rate')
  saveFigures(gcf, fullfile(plot_save_dir, folder_name, 'performance.png'), false)

  % B - frequency change over time
  trial_group_success = cell(1,numel(trial_groups));
  for i=1:numel(trial_group_success)
    [~,trial_group_success{i},~] = getFreqSuccess(trial_groups{i}, 'pow_labs', analysis_idxs, freq_cutoff, params.Fs);
  end
  trial_group_success = cat(3, trial_group_success{:}); % freq x num_splits x num_groups

  y = squeeze(trial_group_success(:,num_splits,:)); % freq x time
  y = y - repmat(suc, 1, size(y,1))'; % calculate delta
  freqs = getFourierFreqs(analysis_idxs, params.Fs);
  freqs = freqs(freqs < freq_cutoff);
  
  figure(plotSpectrogram(y', x, freqs, false))
  colormap jet
  saveFigures(gcf, fullfile(plot_save_dir, folder_name, 'freq_perf_changes.png'), false)
  
  paren = @(A,varargin) A(varargin{:});
  coefs_diff = apply(@(y) paren(corrcoef(conv(y, [1,1]/2, 'valid'), diff(suc)), 2), y, 2, 1);
  figure
  plot(freqs, coefs_diff)
  title('Correlation of harm/benefit of frequency with learning rate')
  saveFigures(gcf, fullfile(plot_save_dir, folder_name, 'freq_corr.png'), false)
end
%% Can a linear discriminant predict trial outcome with just average trace?
% in a word, no
n_folds = 10;

data = single([tbl.trace{:}]); % time x trial
data = [data; (tbl.shape-tbl.noise)'];
[lab, classes] = grp2idx(tbl.result);
disp(classes)
idxs = crossvalind('KFold', size(data,2), n_folds);
n_classes = numel(unique(lab));
confusion_mat = zeros(n_classes, n_classes, n_folds);
eq_group_idxs = cell(1,n_classes);
for i=1:n_folds
  [eq_group_idxs{:}] = equalGroups(lab(idxs~=i));
  train = data(:, idxs~=i);
  train = train(:, vertcat(eq_group_idxs{:}));
  train_labs = lab(idxs~=i);
  train_labs = train_labs(vertcat(eq_group_idxs{:}));
  test = data(:, idxs==i);
  mdl = fitcdiscr(train', train_labs);
  pred_labs = double(categorical(predict(mdl, test')));
  confusion_mat(:,:,i) = accumarray([pred_labs, lab(idxs==i)], 1, [n_classes,n_classes], @sum);
end
disp(mean(confusion_mat,3))
%%
a = spike_tbl.psth(spike_tbl.success,x<0)';
b = spike_tbl.psth(~spike_tbl.success,x<0)';

a1 = mean(fft(a),2);
a1 = a1(1:end/2+1);
a1 = a1(freqs<100);
a2 = fft(mean(a,2));
a2 = a2(1:end/2+1);
a2 = a2(freqs<100);

plot(a1,'.')
hold on
plot(a2, '.')



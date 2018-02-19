% ?? TODO: make phases just be hilbert at shape onset
%% Parameters
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
freq_cutoff = 100;
alpha_freq = 10;
folder_name = 'jan_grant_app';

mt_params = [];
mt_params.tapers = [params.T*params.W, 2*params.T*params.W-1];
mt_params.fpass = [0,params.freq_cutoff];
mt_params.Fs = params.Fs;
mt_params.trialave = true; % this will actually take the mean across all *electrodes* not trials
%% Initialize
% low / high power and success / failure with LFPs and spikes
% HHT on LFPs
% smoothed spikes for each L/H and S/F
fprintf('Selecting data...\n');
spikes = cell(size(conds));
lfps = cell(size(conds));
for i=1:numel(conds)
  spikes{i} = cell(1,numel(dhs));
  lfps{i} = cell(1,numel(dhs));
  for j=1:numel(dhs)
    fprintf('%s, day %d / %d\n', names{i}, j, numel(dhs));
    [data, t] = dhs(j).getDataSlices('spike', select_range, conds{i});
    spikes{i}{j} = data;
    lfps{i}{j} = expandCellArray(dhs(j).electrode_mapping, ...
      dhs(j).getDataSlices('lfp', select_range, ...
       cell2struct({t, false, [0,1]},{'vec','negate','range'},2)));
  end
  valid_days = ~cellfun(@(c) isempty(c) || all(cellfun(@isempty, c)), spikes{i});
  spikes{i} = spikes{i}(valid_days);
  lfps{i} = lfps{i}(valid_days);
end

fprintf('\nDetermining trials with relatively high LFP power...\n');
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
      [psds{i}{j}{k},pre_shape_freqs] = mtspectrumc(trial_lfps{k}(1:-select_range(1),:), mt_params);
    end
    psds{i}{j} = double(cellArray2mat(psds{i}{j}));
  end
end

% get whether a trial had high LFP power
combined_psds = cell(size(psds));
for i=1:numel(combined_psds)
  combined_psds{i} = [psds{i}{:}];
end
is_high_trial = cell(size(conds));
med_lfp_pows = median([combined_psds{:}],2);
for i=1:numel(conds)
  is_high_trial{i} = cellfun(@(c) c > med_lfp_pows, psds{i}, 'UniformOutput', false);
end

fprintf('\nCalculating coherences...\n');
% ignore success / failure conditions
combined_lfps = combineCellArrays('single', lfps{1}{:}, lfps{2}{:});
combined_indices = getCombinedIndices([lfps{:}]);
C = coherencyc(combined_lfps{1}(1:-select_range(1),1), combined_lfps{1}(1:-select_range(1),1), mt_params);
all_cohs = zeros(numel(C), 96, 96);
for i=1:96
  fprintf('%d / %d\n', i, 96);
  all_cohs(:,i,i) = 1;
  for j=i+1:96
    [~,idxs1,idxs2] = intersect(combined_indices{i}, combined_indices{j});
    [C,~,~,~,~,tmp] = coherencyc(combined_lfps{i}(1:-select_range(1),idxs1), combined_lfps{j}(1:-select_range(1),idxs2), mt_params);
    assert(all(tmp == pre_shape_freqs));
    pre_shape_freqs = tmp;
    all_cohs(:,i,j) = C;
    all_cohs(:,j,i) = C;
  end
end
[~,alpha_idx] = min(abs(pre_shape_freqs-alpha_freq));

fprintf('\nCalculating spatial distribution of coherences...\n');
sigmas = getSpatialCohErrors(all_cohs);

fprintf('Calculating phases at shape onset...\n');
phases = cell(1,numel(conds));
for cond_num=1:numel(conds)
  phases{cond_num} = getLfpPhases(lfps{cond_num}, pre_shape_freqs, -select_range(1)+1, [1,-select_range], params.Fs);
end
%% Plot 1 - Success rate vs frequnecy and phase
% A
success_counts = sum([is_high_trial{1}{:}],2);
failure_counts = sum([is_high_trial{2}{:}],2);
success_rate = success_counts ./ (success_counts+failure_counts);
% take binomial (bernoulli) error for success rates
n = size([is_high_trial{1}{:}],2) + size([is_high_trial{2}{:}], 2);
std_err = sqrt(success_rate .* (1-success_rate) / n);

plt = figure('Visible', 'off');
x = pre_shape_freqs(pre_shape_freqs < freq_cutoff);
y = success_rate(pre_shape_freqs < freq_cutoff);
std_err = std_err(pre_shape_freqs < freq_cutoff);
plot(x, y);
hold on
plot(x, [y+std_err, y-std_err], 'r:');
plt.Color = 'white';
xlabel('Frequency (Hz)');
ylabel('Success rate');
title('Success rate by large persence of frequency');
saveFigures(plt, fullfile(plot_save_dir, folder_name, 'pngs/1_a.png'), false);
savefig(plt, fullfile(plot_save_dir, folder_name, 'figs/1_a.fig'));
close(plt);

saveToR(x, y, std_err, data_dir, '1_a');

% B
nbins = 30;
success_counts = histcounts(wrapTo2Pi(angle(phases{1}{alpha_idx})), nbins);
[failure_counts,bins] = histcounts(wrapTo2Pi(angle(phases{2}{alpha_idx})), nbins);
success_rate = success_counts ./ (success_counts+failure_counts);
bins = rad2deg(bins);
% calculate binomial error
n = numel(phases{1}{alpha_idx}) + numel(phases{2}{alpha_idx});
std_err = sqrt(success_rate .* (1-success_rate) / n)';

plt = figure('Visible', 'off');
x = conv(bins, 1/2*[1,1], 'valid')';
y = success_rate';
plot(x, y);
hold on
plot(x, [y+std_err, y-std_err], 'r:');
plt.Color = 'white';
xlabel('Phase Angle');
ylabel('Success rate');
title('Success rate for high alpha 9.8-13.7 Hz');
% title('Success rate for beta 17.6-21.5 Hz');
saveFigures(plt, fullfile(plot_save_dir, folder_name, 'pngs/1_b.png'), false);
savefig(plt, fullfile(plot_save_dir, folder_name, 'figs/1_b.fig'));
close(plt);

saveToR(x, y, std_err, data_dir, '1_b');

% C
success_counts = cellfun(@numel, phases{1});
failure_counts = cellfun(@numel, phases{2});
success_rate = success_counts ./ (success_counts+failure_counts);
% take binomial (bernoulli) error for success rates
n = size([is_high_trial{1}{:}],2) + size([is_high_trial{2}{:}], 2);
std_err = sqrt(success_rate .* (1-success_rate) / n);

plt = figure('Visible', 'off');
x = conv(pre_shape_freqs(pre_shape_freqs < freq_cutoff), 1/2*[1,1], 'valid');
idxs = find(pre_shape_freqs < freq_cutoff);
y = success_rate(idxs(1:end-1));
std_err = std_err(idxs(1:end-1));
plot(x, y);
hold on
plot(x, [y+std_err; y-std_err], 'r:');
plt.Color = 'white';
xlabel('Frequency (Hz)');
ylabel('Success rate');
title('Success rate by presence of IMF with dominant frequency');
saveFigures(plt, fullfile(plot_save_dir, folder_name, 'pngs/1_c.png'), false);
savefig(plt, fullfile(plot_save_dir, folder_name, 'figs/1_c.fig'));
close(plt);

saveToR(x, y, std_err, data_dir, '1_c');
%% Plot 2
% A
breaks = linspace(0, hypot(9,9)*400, 6);
dist_bins = [breaks(1:end-1);breaks(2:end)];
electrode_mapping = cell(96,2);
electrode_mapping(:,1) = num2cell(1:96);

[inter_cohs, num_in_bins, ~, std_err] = calculateInterelecDistCoherence(all_cohs, dist_bins, electrode_mapping);
plt = plotInterelecDistCoherence(inter_cohs, pre_shape_freqs, num_in_bins, dist_bins, 'freq_bin', [0, params.freq_cutoff]);
plt.Color = 'white';
ylim([0,1]);
saveFigures(plt, fullfile(plot_save_dir, folder_name, 'pngs/2_a.png'), false);
savefig(plt, fullfile(plot_save_dir, folder_name, 'figs/2_a.fig'));
close(plt);

saveToR(pre_shape_freqs, inter_cohs, std_err, data_dir, '2_a', 'dist_bins', dist_bins);

% B
plt = figure('Visible', 'off');
plot(pre_shape_freqs, sigmas(:,1));
hold on
plot(pre_shape_freqs, sigmas(:,2:3), 'r:');
plt.Color = 'white';
title('Spatial distribution of coherency');
xlabel('Frequency (Hz)');
ylabel('\sigma (mm)');
saveFigures(plt, fullfile(plot_save_dir, folder_name, 'pngs/2_b.png'), false);
savefig(plt, fullfile(plot_save_dir, folder_name, 'figs/2_b.fig'));
close(plt);

saveToR(pre_shape_freqs, sigmas(:,1), sigmas(:,3)-sigmas(:,1), data_dir, '2_b');
%% Plot 3
% must use pre_shape_freqs here since is_high_tr is determined exclusively
% from pre-shape signal
window_size = 50;
[split_spikes, split_lfps] = splitByLowHigh(spikes, lfps, is_high_trial, alpha_idx);
% A = Success / low       C = Failure / low
% B = Success / high      D = Failure / high
plot_3_names = {'a','b','c','d'};
plot_3_titles = {'Success Low', 'Success High', 'Failure Low', 'Failure High'};
all_smoothed = cell(size(split_spikes));
for i=1:numel(split_spikes)
  all_smoothed{i} = cell(size(split_spikes{i}));
  for j=1:numel(split_spikes{i}) % days
    % if a day doesn't have an electrode, smoothed{n} is just set to [], which
    % is later ignored by cellArray2mat(all_smoothed{i})
    [smoothed,t] = cellfun(@(c) calculateSmoothedSpikes(c,gausswin(window_size),select_range(1),diff(select_range)),...
                 split_spikes{i}{j}, 'UniformOutput', false);
    % smoothed is units -> time x trials
    smoothed = reshape([smoothed{:}], [size(smoothed{1}), numel(smoothed)]); % time x trials x units
    all_smoothed{i}{j} = mean(smoothed,3); % time x trials
  end
  all_smoothed{i} = [all_smoothed{i}{:}];
  [plt,x,y,std_err] = plotMeanAndStds(all_smoothed{i} * params.Fs / window_size, 'x', t{1});
  plt.Color = 'white';
  title(sprintf('%s Power', plot_3_titles{i}));
  xlabel('Time (ms)');
  ylabel('Mean firing rate (Hz)');
  saveFigures(plt, fullfile(plot_save_dir, folder_name, sprintf('pngs/3_%s.png', plot_3_names{i})), false);
  savefig(plt, fullfile(plot_save_dir, folder_name, sprintf('figs/3_%s.fig', plot_3_names{i})));
  close(plt);
  
  saveToR(x, y, std_err, data_dir, sprintf('3_%s', plot_3_names{i}));
end
%% Plot 4
plot_4_mt_params = mt_params;
plot_4_mt_params.err = [1, 0.95];
for i=1:numel(all_smoothed)
  % TODO: clean up x<0   code
  [S,f,Serr] = mtspectrumc(all_smoothed{1}(x<0,:), plot_4_mt_params);
  plt = figure('Visible', 'off');
  plot(f,10*log10(S));
  hold on;
  plot(f, 10*log10(Serr), 'r:');
  plt.Color = 'white';
  xlabel('Frequency (Hz)');
  ylabel('Power');
  saveFigures(plt, fullfile(plot_save_dir, folder_name, sprintf('pngs/4_%s.png', plot_3_names{i})), false);
  savefig(plt, fullfile(plot_save_dir, folder_name, sprintf('figs/4_%s.fig', plot_3_names{i})));
  close(plt);
  
  saveToR(f, S, Serr, data_dir, sprintf('4_%s', plot_3_names{i}));
end
%%
kl_divs = zeros(1,numel(pre_shape_freqs)-1);
for i=1:numel(pre_shape_freqs)-1
  kl_divs(i) = KLDiv(histcounts(wrapTo2Pi(angle(phases{1}{i})), 12), ...
                     histcounts(wrapTo2Pi(angle(phases{2}{i})), 12));
end
%%
% total_n = zeros(size(split_spikes));
% for i=[1,2]
%   for j=[1,2]
%     total_n(j,i) = sum(cellfun(@(c) numel(c{1}), split_spikes{j,i}));
%   end
% end

total_n = [];
total_n.low_success = sum(cellfun(@(c) numel(c{1}), split_spikes{1,1}));
total_n.high_success = sum(cellfun(@(c) numel(c{1}), split_spikes{2,1}));
total_n.low_failure = sum(cellfun(@(c) numel(c{1}), split_spikes{1,2}));
total_n.high_failure = sum(cellfun(@(c) numel(c{1}), split_spikes{2,2}));
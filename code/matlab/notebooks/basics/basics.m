monkey_name = 'jaws';
analysis_type = 'preshape';
folder_name = 'basics';

freq_cutoff = 500;

% sorting into power bins
num_splits = 2;

% moving average params
ma_window_size = 100;
ma_step_size = 50; 

% coherence distance params
dist_edges = linspace(0, hypot(9,9)*400, 6);

if strcmp(analysis_type, 'preshape')
  analysis_idxs = 1:256;
  trial_start = -256;
elseif strcmp(analysis_type, 'prenoise')
  analysis_idxs = 1:500;
  trial_start = -500;
end

trial_dir = fullfile(data_dir, 'trials', monkey_name);
lfp_dir = fullfile(data_dir, 'lfps', monkey_name);
dh_dir = fullfile(cache_dir, 'data_handlers', monkey_name);
save_dir = fullfile(cache_dir, 'trial_tables', analysis_type, monkey_name);
all_save_dir = fullfile(cache_dir, 'trial_tables', 'all', monkey_name);

mt_params = [];
mt_params.tapers = dpss(numel(analysis_idxs), params.T*params.W, 2*params.T*params.W-1);
mt_params.Fs = params.Fs;
mt_params.trialave = true;
mt_params.pad = -1; % don't pad

if strcmp(monkey_name, 'jaws')
  dates = datetime('2013-05-1'):datetime('2013-07-19');
  learn_start = datetime('2013-06-3');
else
  dates = datetime('2011-7-22'):datetime('2011-12-20');
  learn_start = datetime('2011-11-7');
end
%%
all_tbl = makeTrialTable(dates, trial_dir, lfp_dir, dh_dir, save_dir, analysis_type, ...
  'clean', true, ...
  'min_reaction_time', params.length_min_reaction_time, ...
  'monkey_name', monkey_name);
%%
all_tbl.success = strcmp(all_tbl.result, 'true_positive');
all_tbl.trace = cellfun(@(c) mean(c,2), all_tbl.lfps, 'UniformOutput', false);
all_tbl.psd = getPsd(all_tbl, analysis_idxs, mt_params);

grouped = groupBy(all_tbl, 'date');
cons = cell(1, numel(grouped));
parfor ii=1:numel(grouped)
  fprintf("Connectivity: day %s\n", grouped{ii}(1,:).date{1});
  % connecitivy uses PPC
  % con = chan x chan x freq
  con = connectivityLFP(grouped{ii}.lfps, trial_start, grouped{ii}.shape);
  cons{ii} = mean(reshape(con, [size(con,1)*size(con,2), size(con,3)]),1);
end
cons = vertcat(cons{:})'; % freq x days

% NOTE: false alarms (false positives) don't make sense for shape-aligned
% data -- only makes sense with a prenoise dataset
tbl = all_tbl(contains(all_tbl.result, {'true_positive', 'false_negative', 'false_positive'}),:);
%% Learning Rates
% sensitivity and false discovery rate
tp = tbl.success;
fn = strcmp(tbl.result, 'false_negative');
fp = strcmp(tbl.result, 'false_positive');
rates = movmean([tp, fn, fp], 300, 1, 'EndPoints', 'discard');
rates = rates(1:ma_step_size:end,:);
x = movmean(1:size(tbl,1), 300, 'EndPoints', 'discard');
x = x(1:ma_step_size:end);
figure
h = plot(x', rates);
hold on
learn_idx = x(find(datetime(tbl.date(round(x))) >= learn_start, 1));
plot([learn_idx,learn_idx], [0, 1], 'k--', 'LineWidth', 1)
legend(h, {'Success', 'Detection failure', 'False alarm'}, 'Location', 'northoutside')
saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('learning_rates_%s.png', monkey_name)));
%% Traces
figure
lfps = [tbl.trace{:}]; % time x trial
lfps = movmean(lfps, ma_window_size, 2, 'EndPoints', 'discard');
lfps = lfps(:,1:ma_step_size:end);

imagesc(flipud(lfps))
xlabel('Trial Index')
ylabel('Time relative to shape (ms)');
yticks(1:50:size(lfps,1));
yticklabels(arrayfun(@num2str, analysis_idxs(end)+trial_start:-50:trial_start+1, 'UniformOutput', false))
colormap jet
colorbar
x = movmean(1:size(tbl,1), ma_window_size, 'EndPoints', 'discard');
t = x(1:ma_step_size:end);
times = linspace(t(1),t(end),10);
time_ticks = interp1(times, linspace(1,numel(t),10), times);
time_labels = arrayfun(@(x) sprintf('%.0f',x), times, 'UniformOutput', false);
set(gca, 'XTick', time_ticks, 'XTickLabel', time_labels);
xlabel('Trail Index');
hold on
learn_idx = find(datetime(tbl.date(round(t))) == learn_start, 1);
plot([learn_idx,learn_idx], [1,numel(analysis_idxs)], 'k--', 'LineWidth', 3)
saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('all_traces_shape_%s.png', monkey_name)));

lfps = [tbl.trace{:}];
learn_idx = find(datetime(tbl.date) >= learn_start, 1);
idxs1 = round(linspace(1, learn_idx-1, 5));
idxs2 = round(linspace(learn_idx, size(tbl,1), 5));

lim = [min(lfps(:)), max(lfps(:))] * 1.05;
figure
plot(analysis_idxs'+trial_start, lfps(:,idxs1))
title('Before learning')
xlabel('Time (ms)')
ylabel('Potential (mV)')
ylim(lim)
saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('sample_traces_before_%s.png', monkey_name)));

figure
plot(analysis_idxs'+trial_start, lfps(:,idxs2))
title('After learning')
xlabel('Time (ms)')
ylabel('Potential (mV)')
ylim(lim)
saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('sample_traces_after_%s.png', monkey_name)));

idxs=analysis_idxs+trial_start>=0;
figure(plotMeanAndStds(lfps(idxs,(1:end) < learn_idx), 'x', analysis_idxs(idxs)'+trial_start))
title('Before learning')
xlabel('Time (ms)')
ylabel('Potential (mV)')
saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('ave_traces_before_%s.png', monkey_name)));

figure(plotMeanAndStds(lfps(idxs,(1:end) >= learn_idx), 'x', analysis_idxs(idxs)'+trial_start))
title('After learning')
xlabel('Time (ms)')
ylabel('Potential (mV)')
saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('ave_traces_after_%s.png', monkey_name)));

% sample trials
idxs = [1000, 1031, 2000, 4000, 5000, 6000];
for i=1:numel(idxs)
  figure
  plot(analysis_idxs'+trial_start, tbl.lfps{idxs(i)});
  xlabel('Time (ms)')
  ylabel('Potential (mV)')
  title(sprintf('Trial #%d', idxs(i)))
  saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('sample_trials_%d_%s.png', i, monkey_name)));
end
%% PSD
x = movmean(1:size(tbl,1), ma_window_size, 'EndPoints', 'discard');
x = x(1:ma_step_size:end);
psd = movmean(tbl.psd, ma_window_size, 1, 'EndPoints', 'discard'); % trial x freq
psd = psd(1:ma_step_size:end,:);

figure(plotSpectrogram(psd, x, getFourierFreqs(analysis_idxs,1000)));
colormap jet
times = linspace(x(1),x(end),10);
time_ticks = interp1(times, linspace(1,numel(x),10), times);
time_labels = arrayfun(@(x) sprintf('%.0f',x), times, 'UniformOutput', false);
set(gca, 'XTick', time_ticks, 'XTickLabel', time_labels);
xlabel('Trail Index');
hold on
learn_idx = find(datetime(tbl.date(round(x))) == learn_start, 1);
plot([learn_idx,learn_idx], [1,numel(analysis_idxs)], 'k--', 'LineWidth', 1.5)
saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('spectra_%s.png', monkey_name)));
%% Phase locking
% for now (speed / memory reasons) this is only by day, not trial groupings
x = 1:size(cons,2);

figure(plotSpectrogram(cons', x, getFourierFreqs(analysis_idxs,1000), false));
colormap jet
times = linspace(x(1),x(end),10);
time_ticks = interp1(times, linspace(1,numel(x),10), times);
time_labels = arrayfun(@(x) sprintf('%.0f',x), times, 'UniformOutput', false);
set(gca, 'XTick', time_ticks, 'XTickLabel', time_labels);
xlabel('Day Number');
hold on
dates = sort(datetime(unique(tbl.date)));
learn_idx = find(dates == learn_start, 1);
plot([learn_idx,learn_idx], [1,numel(analysis_idxs)], 'k--', 'LineWidth', 1.5)
saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('phase_locking_%s.png', monkey_name)));
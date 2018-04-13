monkey_name = 'jaws';
analysis_type = 'noise';
folder_name = 'capstone';

alpha_range = [7.5,25];
gamma_range = [31,100];
all_ranges = [alpha_range; gamma_range];

num_splits = 3;

trial_dir = fullfile(data_dir, 'trials', monkey_name);
lfp_dir = fullfile(data_dir, 'lfps', monkey_name);
dh_dir = fullfile(cache_dir, 'data_handlers', monkey_name);
save_dir = fullfile(cache_dir, 'trial_tables', analysis_type, monkey_name);
output_dir = fullfile(cache_dir, 'output', monkey_name);

mt_params = [];
mt_params.tapers = [params.T*params.W, 2*params.T*params.W-1];
mt_params.Fs = params.Fs;
mt_params.trialave = true;

if strcmp(monkey_name, 'jaws')
  dates = datetime('2013-06-3'):datetime('2013-07-19');
else
  dates = datetime('2011-11-7'):datetime('2011-12-20');
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
tbl = tbl(tbl.shapecoh < 70,:); % hardest 2
grouped = groupBy(tbl, 'date');
%% Coherence by day
% day_cohs = cell(1,numel(grouped));
% coh_distribution = cell(1,numel(grouped));
% parfor ii=1:numel(grouped)
%   fprintf("Day %s\n", grouped{ii}(1,:).date{1});
%   day_cohs{ii} = cohByDay(grouped{ii}.lfps, grouped{ii}(1,:).electrode_mapping{1}, analysis_idxs, mt_params);
%   %coh_distribution{ii} = getSpatialCohErrors(day_cohs{ii});
% end
% %% Saving coherence by day
% x = unique(tbl.date);
% freqs = getFourierFreqs(analysis_idxs,1000);
% y = cat(3, coh_distribution{:});
% y = squeeze(y(:,1,:)); % just get the mean sigmas, not the std errs
% saveToR(x, y, [], data_dir, folder_name, sprintf('coh_%s_day_%s', analysis_type, monkey_name), 'freqs', freqs);
%% LFP phases
if strcmp(analysis_type, 'shape')
  pows = rowfun(@(x) {mtspectrumc(mean(x{1},2), mt_params)}, tbl, 'InputVariables', {'lfps'});
  pows = rowfun(@(x) {x{1}/max(x{1})}, pows);
  pows = [pows.Var1{:}];
  binned_pows = zeros(size(all_ranges,1),size(pows,2));
  freqs = getFourierFreqs(analysis_idxs,1000);
  for i=1:size(all_ranges,1)
    idxs = between(freqs, all_ranges(i,1), all_ranges(i,2));
    binned_pows(i,:) = mean(pows(idxs,:),1);
  end
  pow_bins = quantile(binned_pows, (0:num_splits)/num_splits, 2);
  pow_labs = apply(@(i) discretize(binned_pows(i,:), pow_bins(i,:)), ...
                   1:size(all_ranges,1), 1, size(binned_pows,2)); % discretize each row (row = trial)
  tbl.pow_labs = pow_labs;
  
  phases = rowfun(@(x) {angle(fft(mean(x{1},2)))}, tbl, 'InputVariables', {'lfps'});
  phases = rowfun(@(x) {x{1}(1:floor(end/2)+1)}, phases);
  tbl.phases = phases.Var1;
  [success_groups,groups] = groupBy(tbl, 'result');
  success_idx = find(strcmp(groups, 'true_positive'));
  failure_idx = find(strcmp(groups, 'false_negative'));
  perf_idxs = [success_idx,failure_idx];
  grouped_phases = cell(2, size(all_ranges,1));
  for i=1:size(all_ranges,1)
    idxs = between(freqs, all_ranges(i,1), all_ranges(i,2));
    for j=1:numel(perf_idxs)
      high_pow_tbl = success_groups{perf_idxs(j)};
      high_pow_tbl = high_pow_tbl(high_pow_tbl.pow_labs(:,i)==num_splits,:);
      grouped_phases{j,i} = rowfun(@(x) {x{1}(idxs)}, high_pow_tbl, 'InputVariables', {'phases'});
      grouped_phases{j,i} = vertcat(grouped_phases{j,i}.Var1{:});
    end
  end
  saveToR(grouped_phases, [], [], data_dir, folder_name, sprintf('bigger_bins_phases_%s', monkey_name));
end
%% Reaction time vs power bin for low and high frequencies
if strcmp(analysis_type, 'shape')
  s_tbl = tbl(strcmp(tbl.result, 'true_positive'),:);
  rxn_times = s_tbl.saccade - s_tbl.shape;
  % sort into groups based on pow_lab
  rxn_groups_low = accumarray(s_tbl.pow_labs(:,1), rxn_times, [], @(x) {x});
  rxn_groups_high = accumarray(s_tbl.pow_labs(:,2), rxn_times, [], @(x) {x});
  saveToR(1:num_splits, {rxn_groups_low, rxn_groups_high}, [], ...
          data_dir, folder_name, sprintf('rxn_time_pow_bin_%s', monkey_name));
  saveToR([], [], [], data_dir, folder_name, sprintf('1_d_%s', monkey_name), 'rxn_time', rxn_times);
end
%% Coherence by trial
% grouping by day necessary for doing coherence with one consistent electrode mapping
breaks = linspace(0, hypot(9,9)*400, 6);
freqs = getFourierFreqs(analysis_idxs,1000);
[~,freq_idx] = min(abs(freqs - 20));
window_size = 300;
step_size = 50;

trial_cohs = cell(1,numel(grouped));
trial_mt_params = mt_params;
trial_mt_params.trialave = false;
parfor ii=1:numel(grouped)
  fprintf("Day %s\n", grouped{ii}(1,:).date{1});
  trial_cohs{ii} = cohByDay(grouped{ii}.lfps, grouped{ii}(1,:).electrode_mapping{1}, analysis_idxs, trial_mt_params, freq_idx);
end
trial_cohs = cat(4, trial_cohs{:});
[inter_cohs,~,~,std_err] = groupByDist(trial_cohs, breaks);
farthest = squeeze(inter_cohs(:,:,end)); % freq x trials
std_err = squeeze(std_err(:,:,end));

% TODO: perform better moving window calculation that takes operates on
% binned_cohs
x = movmean(1:size(farthest,2), window_size, 'EndPoints', 'discard');
x = x(1:step_size:end);
y = movmean(farthest, window_size, 'EndPoints', 'discard');
y = y(1:step_size:end);
std_err = movmean(std_err, window_size, 'EndPoints', 'discard');
std_err = std_err(1:step_size:end);
saveToR(x', y', std_err', data_dir, folder_name, sprintf('trial_coh_%s_%s', analysis_type, monkey_name));

% Sucess rate
success = strcmp(tbl.result, 'true_positive');
x = movmean(1:size(tbl,1), window_size, 'EndPoints', 'discard');
x = x(1:step_size:end);
y = movmean(success', window_size, 'EndPoints', 'discard');
y = y(1:step_size:end);
saveToR(x', y', [], data_dir, folder_name, sprintf('success_rate_coh_%s_%s', analysis_type, monkey_name));
%%
% ds = datastore(save_dir, 'Type', 'file', 'FileExtensions', '.mat', ...
%   'ReadFcn', @(filename) getfield(load(filename), 'tbl'));
% 
% mr = mapreducer;
% outds = mapreduce(ds, @(x,y,z) analysisByRow(x,y,z,params), @reduceByDay, mr, ...
%   'OutputFolder', output_dir);
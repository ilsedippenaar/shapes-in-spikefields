%% Initialize
names = {'success', 'failure'};
pow_names = {'low', 'high'};
freq_range = [17,22];
select_range = [-256, 256];
folder_name = sprintf('%d_to_%d', freq_range(1), freq_range(2));

success_shape_cond = cell2struct([...
  {'shape', false, [0,1]}; ...
  {'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]}; ...
  {'saccade', true, [-1000, params.length_min_reaction_time]}; ... % shapes can come on *after* a saccade too
  {'saccade', false, [params.length_min_reaction_time, 1000]}], {'vec','negate','range'}, 2);

failure_shape_cond = cell2struct([...
  {'shape', false, [0,1]}; ...
  {'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]}; ...
  {'saccade', true, [-1000, 1000]}], {'vec','negate','range'}, 2);

res_conds = {success_shape_cond, failure_shape_cond}; % res = trial result
res_spikes = cell(size(res_conds));
res_times = cell(size(res_conds));
res_lfps = cell(size(res_conds));
combined_indices = cell(size(res_conds));
for i=1:numel(res_conds)
  res_spikes{i} = cell(1,numel(dhs));
  res_times{i} = cell(1,numel(dhs));
  res_lfps{i} = cell(1,numel(dhs));
  num_in_day = zeros(1,numel(dhs));
  for j=1:numel(dhs)
     [data, t] = dhs(j).getDataSlices('spike', select_range, res_conds{i});
     res_spikes{i}{j} = data;
     res_times{i}{j} = t;
     res_lfps{i}{j} = expandCellArray(dhs(j).electrode_mapping, ...
       dhs(j).getDataSlices('lfp', select_range, ...
        cell2struct({t, false, [0,1]},{'vec','negate','range'},2)));
      num_in_day(j) = numel(data{1}); % number of trials
  end
  valid_days = ~cellfun(@(c) all(cellfun(@isempty, c)), res_spikes{i});
  res_spikes{i} = res_spikes{i}(valid_days);
  res_lfps{i} = res_lfps{i}(valid_days);
  
  % TODO: this is super wrong - probably fixed now?
  num_in_day = num_in_day(valid_days);
  tmp = cumsum(num_in_day);
  trial_idxs_bounds = [1,tmp(1:end-1)+1; tmp];
  all_idxs = cell(size(trial_idxs_bounds,2));
  for j=1:size(trial_idxs_bounds,2)
    all_idxs{j} = trial_idxs_bounds(1,j):trial_idxs_bounds(2,j);
  end
  combined_indices{i} = cell(1,96);
  for j=1:96
    combined_indices{i}{j} = cellArray2mat(all_idxs(cellfun(@(c) ~isempty(c{j}), res_lfps{i})));
  end
end

res_mean_pows = cell(size(res_lfps));
% calculate average PSD for each trial
mt_params = [];
mt_params.tapers = [params.T*params.W, 2*params.T*params.W-1];
mt_params.Fs = 1000;
mt_params.trialave = true;
res_psd = cell(size(res_lfps));
for i=1:numel(res_lfps) % result
  res_psd{i} = cell(size(res_lfps{i}));
  for j=1:numel(res_lfps{i}) % day
    fprintf('%d\n', j);
    % rearrange to have electrodes nested in a trial
    tmp = cell(1,numel(res_lfps{i}{j}));
    for k=1:numel(res_lfps{i}{j}) % electrode
      % trials to cell array
      tmp{k} = mat2cell(res_lfps{i}{j}{k}, size(res_lfps{i}{j}{k},1), repelem(1,size(res_lfps{i}{j}{k},2)));
    end
    trial_lfps = combineCellArrays('single', tmp{:});
    res_psd{i}{j} = cell(1,numel(trial_lfps));
    for k=1:numel(trial_lfps)
      [res_psd{i}{j}{k},f] = mtspectrumc(trial_lfps{k}(1:-select_range(1),:), mt_params);
    end
    res_psd{i}{j} = double(cellArray2mat(res_psd{i}{j}));
  end
  % get whether a trial had high LFP power
  freq_idxs = and(f >= freq_range(1), f < freq_range(2));
  res_mean_pows{i} = cellArray2mat(cellfun(@(c) mean(c(freq_idxs,:)), res_psd{i}, 'UniformOutput', false));
end
is_high_tr = cell(size(res_conds));
med_lfp_pow = median(cellArray2mat(res_mean_pows));
for i=1:numel(res_lfps)
  is_high_tr{i} = cellfun(@(c) mean(c(freq_idxs,:),1) > med_lfp_pow, res_psd{i}, 'UniformOutput', false);
end

% divide spikes and lfps into four groups
res_pow_spikes = cell(2,2); % S/F by low/high power
res_pow_lfps = cell(size(res_pow_spikes));
lfp_pow_names = cell(size(res_pow_spikes));
for i=[1,2] % S/F
  for j=[1,2] % LFP power
    res_pow_spikes{i,j} = cell(size(res_spikes{i}));
    res_pow_lfps{i,j} = cell(size(res_spikes{i}));
    lfp_pow_names{i,j} = sprintf('%s_%s_power', names{i}, pow_names{j});
    for k=1:numel(res_spikes{i}) % day
      res_pow_spikes{i,j}{k} = cell(size(res_spikes{i}{k}));
      for l=1:numel(res_spikes{i}{k}) % unit
        res_pow_spikes{i,j}{k}{l} = res_spikes{i}{k}{l}(~xor(j-1, is_high_tr{i}{k})); % 0 negates is_high_tr, 1 leaves it the same
      end
      res_pow_lfps{i,j}{k} = cell(size(res_lfps{i}{k}));
      for l=1:numel(res_lfps{i}{k}) % electrode
        if ~isempty(res_lfps{i}{k}{l})
          res_pow_lfps{i,j}{k}{l} = res_lfps{i}{k}{l}(:,~xor(j-1, is_high_tr{i}{k}));
        end
      end
    end
  end
end
%% Plot PSTHs split by S/F and LFP power
all_smoothed = cell(size(res_pow_spikes));
for i=1:numel(res_pow_spikes)
  all_smoothed{i} = cell(size(res_pow_spikes{i}));
  for j=1:numel(res_pow_spikes{i}) % days
    % if a day doesn't have an electrode, tmp{n} is just set to [], which
    % is later ignored by cellArray2mat(all_smoothed{i})
    tmp = cellfun(@(c) calculateSmoothedSpikes(c,gausswin(30),select_range(1),diff(select_range)),...
                 res_pow_spikes{i}{j}, 'UniformOutput', false);
    all_smoothed{i}{j} = cellArray2mat(tmp);
  end
  all_smoothed{i} = cellArray2mat(all_smoothed{i});
  plt = plotMeanAndStds(all_smoothed{i},'x', select_range(1):select_range(2)-1);
  saveFigures(plt, fullfile(plot_save_dir, 'spikes', folder_name, sprintf('%s.png', lfp_pow_names{i})), false);
  savefig(plt, fullfile(plot_save_dir, 'figs/spike_field', folder_name, sprintf('%s.fig', lfp_pow_names{i})));
  close(plt);
end
%%
breaks = linspace(0, hypot(9,9)*400, 6);
dist_bins = [breaks(1:end-1);breaks(2:end)];
electrode_mapping = cell(96,2);
electrode_mapping(:,1) = num2cell(1:96);

mt_params = [];
mt_params.tapers = [params.T*params.W, 2*params.T*params.W-1];
mt_params.fpass = [0,params.freq_cutoff];
mt_params.Fs = params.Fs;
mt_params.trialave = true;
all_cohs = cell(1,numel(res_conds));
for cond_num=1:numel(res_conds)
  combined_lfps = combineCellArrays('single', res_lfps{cond_num}{:});
  C = coherencyc(combined_lfps{j}(:,1), combined_lfps{j}(:,1), mt_params);
  all_cohs{cond_num} = zeros(numel(C), 96, 96);
  for j=1:96
    fprintf('%d / %d\n', j, 96);
    all_cohs{cond_num}(:,j,j) = 1;
    for k=j+1:96
      [~,idxs1,idxs2] = intersect(combined_indices{cond_num}{j},combined_indices{cond_num}{k});
      [C,~,~,~,~,f] = coherencyc(combined_lfps{j}(:,idxs1), combined_lfps{k}(:,idxs2), mt_params);
      all_cohs{cond_num}(:,j,k) = C;
      all_cohs{cond_num}(:,k,j) = C;
    end
  end
end
%% Plot intereclectrode distance coherence
for cond_num=1:numel(res_conds)
  [inter_cohs, num_in_bins, binned_cohs] = calculateInterelecDistCoherence(all_cohs{cond_num}, dist_bins, electrode_mapping);
  plt = plotInterelecDistCoherence(inter_cohs, f, num_in_bins, dist_bins, 'freq_bin', [0, params.freq_cutoff]);
  plt.Color = 'white';
  saveFigures(plt, fullfile(plot_save_dir, 'coherence', folder_name, sprintf('interelectrode_dist_%s.png', names{cond_num})), false);
  savefig(plt, fullfile(plot_save_dir, 'figs/coherence', folder_name, sprintf('%s.fig', names{cond_num})));
  close(plt);
end
%% Count number in each group
total_n = zeros(size(res_pow_spikes));
for i=[1,2]
  for j=[1,2]
    total_n(i,j) = sum(cellfun(@(c) numel(c{1}), res_pow_spikes{i,j}));
  end
end
%% Initialize
names = {'tp', 'fn'};
freq_range = [12,15];
select_range = [-250, 250];

tp_shape_cond = cell2struct([...
  {'shape', false, [0,1]}; ...
  {'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]};...
  {'saccade', false, [0, 1000]}], {'vec','negate','range'}, 2);

fn_shape_cond = cell2struct([...
  {'shape', false, [0,1]}; ...
  {'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]};...
  {'saccade', true, [0, 1000]}], {'vec','negate','range'}, 2);

res_conds = {tp_shape_cond, fn_shape_cond}; % res = trial result
res_spikes = cell(size(res_conds));
res_times = cell(size(res_conds));
res_lfps = cell(size(res_conds));
for i=1:numel(res_conds)
  res_spikes{i} = cell(1,numel(dhs));
  res_times{i} = cell(1,numel(dhs));
  res_lfps{i} = cell(1,numel(dhs));
  for j=1:numel(dhs)
     [data, t] = dhs(j).getDataSlices('spike', select_range, res_conds{i});
     res_spikes{i}{j} = data;
     res_times{i}{j} = t;
     res_lfps{i}{j} = dhs(j).getDataSlices('lfp', select_range, ...
       cell2struct({t, false, [0,1]},{'vec','negate','range'},2));
  end
  valid_days = ~cellfun(@(c) all(cellfun(@isempty, c)), res_spikes{i});
  res_spikes{i} = res_spikes{i}(valid_days);
  res_lfps{i} = res_lfps{i}(valid_days);
end

is_high_tr = cell(size(res_conds));
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
  % get power in frequency range and use it to divide spikes into two groups
  freq_idxs = and(f >= freq_range(1), f < freq_range(2));
  med_lfp_pow = median(cellArray2mat(cellfun(@(c) mean(c(freq_idxs,:)), res_psd{i}, 'UniformOutput', false)));
  is_high_tr{i} = cellfun(@(c) mean(c(freq_idxs,:),1) > med_lfp_pow, res_psd{i}, 'UniformOutput', false);
end
%%
res_pow_spikes = cell(2,2);
for i=[1,2] % S/F
  for j=[1,2] % LFP power
    res_pow_spikes{i,j} = cell(size(res_spikes{i}));
    for k=1:numel(res_spikes{i}) % day
      res_pow_spikes{i,j}{k} = cell(size(res_spikes{i}{k}));
      for l=1:numel(res_spikes{i}{k}) % unit
        res_pow_spikes{i,j}{k}{l} = res_spikes{i}{k}{l}(~xor(j-1, is_high_tr{i}{k})); % 0 negates is_high_tr, 1 leaves it the same
      end
    end
  end
end

% plot (smoothed) PSTH
lfp_pow_names = {...
  sprintf('success_low_power_%d_to_%d', freq_range(1), freq_range(2)),...
  sprintf('success_high_power_%d_to_%d', freq_range(1), freq_range(2));
  sprintf('failure_low_power_%d_to_%d', freq_range(1), freq_range(2)),...
  sprintf('failure_high_power_%d_to_%d', freq_range(1), freq_range(2))};
all_smoothed = cell(size(res_pow_spikes));
for i=1:numel(res_pow_spikes)
  all_smoothed{i} = cell(size(res_pow_spikes{i}));
  for j=1:numel(res_pow_spikes{i})
    tmp = cellfun(@(c) calculateSmoothedSpikes(c,gausswin(50),select_range(1),diff(select_range)),...
                 res_pow_spikes{i}{j}, 'UniformOutput', false);
    all_smoothed{i}{j} = cellArray2mat(tmp);
  end
  all_smoothed{i} = cellArray2mat(all_smoothed{i});
  plt = plotMeanAndStds(all_smoothed{i},'x', select_range(1):select_range(2)-1);
  saveFigures(plt, fullfile(plot_save_dir, 'all_days/spikes', sprintf('psth_shape_%s.png', lfp_pow_names{i})));
end
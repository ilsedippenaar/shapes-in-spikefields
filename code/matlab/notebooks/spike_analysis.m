% PSTH across H
% PSTH across H by spectral LFP power in a frequency
% ISI vs time
%% Initialize
select_range = [-250, 250];
shape_cond = cell2struct([...
  {'shape', false, [0,1]}; ...
  {'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]}], {'vec','negate','range'}, 2);

spikes = cell(1,numel(dhs));
lfps = cell(1,numel(dhs));
for i=1:numel(dhs)
  fprintf('%d\n',i);
  spikes{i} = dhs(i).getDataSlices('spike', select_range, shape_cond);
  lfps{i} = dhs(i).getDataSlices('lfp', select_range, shape_cond);
end
% exclude days without any data and days where, due to recording, there
% are more trials with lfps than with spikes, just to be safe
valid_days = and(~cellfun(@(c) all(cellfun(@isempty, c)), spikes), ...
  cellfun(@(c) size(c{1},2), spikes) == cellfun(@(c) size(c{1},2), lfps));
spikes = spikes(valid_days);
lfps = lfps(valid_days);
momnames = {'tp', 'fn'};
tp_shape_cond = cell2struct([...
  {'shape', false, [0,1]}; ...
  {'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]};...
  {'saccade', false, [0, 1000]}], {'vec','negate','range'}, 2);

fn_shape_cond = cell2struct([...
  {'shape', false, [0,1]}; ...
  {'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]};...
  {'saccade', true, [0, 1000]}], {'vec','negate','range'}, 2);

tr_shape_conds = {tp_shape_cond, fn_shape_cond}; % tr = trial result
tr_spikes = cell(size(tr_shape_conds));
for i=1:numel(tr_shape_conds)
  tr_spikes{i} = cell(1,numel(dhs));
  for j=1:numel(dhs)
    tr_spikes{i}{j} = dhs(j).getDataSlices('spike', select_range, tr_shape_conds{i});
  end
  tr_spikes{i} = tr_spikes{i}(and(valid_days, ...
    ~cellfun(@(c) all(cellfun(@isempty, c)), tr_spikes{i})));
end

%% LFP power initialize
% calculate average PSD for each trial
mt_params = [];
mt_params.tapers = [params.T*params.W, 2*params.T*params.W-1];
mt_params.Fs = 1000;
mt_params.trialave = true;
trial_psd = cell(size(lfps));
for i=1:numel(lfps)
  tmp = cell(1,numel(lfps{i}));
  for j=1:numel(lfps{i})
    tmp{j} = mat2cell(lfps{i}{j}, size(lfps{i}{j},1), repelem(1,size(lfps{i}{j},2)));
  end
  trial_lfps = combineCellArrays('single', tmp{:});
  trial_psd{i} = cell(1,numel(trial_lfps));
  for j=1:numel(trial_lfps)
    trial_psd{i}{j} = mtspectrumc(trial_lfps{j}, mt_params);
  end
  trial_psd{i} = double(cellArray2mat(trial_psd{i}));
end

% get power in frequency range and use it to divide spikes into two groups
freq_range = [12,15];
freq_idxs = and(f >= freq_range(1), f < freq_range(2));
med_lfp_pow = median(cellArray2mat(cellfun(@(c) mean(c(freq_idxs,:)),trial_psd, 'UniformOutput', false)));
is_high_tr = cellfun(@(c) mean(c(freq_idxs,:),1) > med_lfp_pow, trial_psd, 'UniformOutput', false);
%% PSTH (all neurons)
smoothed = cell(size(spikes));
for i=1:numel(spikes)
  smoothed{i} = cellArray2mat(...
    cellfun(@(c) calculateSmoothedSpikes(...
                  c,gausswin(50),select_range(1),diff(select_range)),...
      spikes{i}, 'UniformOutput', false));
end
plt = plotMeanAndStds(cellArray2mat(smoothed),'x', select_range(1):select_range(2)-1);
saveFigures(plt, fullfile(plot_save_dir, 'all_days/spikes', 'psth_shape_all.png'));
%% PSTH by succcess/failure (TP vs FN)
for i=1:numel(tr_spikes)
  tr_smoothed_spikes = cell(size(tr_spikes{i}));
  for j=1:numel(tr_spikes{i})
    tmp = cellfun(@(c) calculateSmoothedSpikes(c,gausswin(50),select_range(1),diff(select_range)),...
                 tr_spikes{i}{j}, 'UniformOutput', false);
    tr_smoothed_spikes{j} = cellArray2mat(tmp);
  end
  plt = plotMeanAndStds(cellArray2mat(tr_smoothed_spikes),'x', select_range(1):select_range(2)-1);
  saveFigures(plt, fullfile(plot_save_dir, 'all_days/spikes', sprintf('psth_shape_%s.png', names{i})));
end
%% PSTH by LFP power
% divide the spikes
lfp_pow_spikes = cell(1,2);
for i=[1,2] % low lfp power then high
  lfp_pow_spikes{i} = cell(size(spikes));
  for j=1:numel(spikes)
    lfp_pow_spikes{i}{j} = cell(size(spikes{j}));
    for k=1:numel(spikes{j})
      lfp_pow_spikes{i}{j}{k} = spikes{j}{k}(~xor(i-1, is_high_tr{j})); % 0 negates is_high_tr, 1 leaves it the same
    end
  end
end

% plot (smoothed) PSTH
lfp_pow_names = {sprintf('low_power_%d_to_%d', freq_range(1), freq_range(2)),...
  sprintf('high_power_%d_to_%d', freq_range(1), freq_range(2))};
for i=1:numel(lfp_pow_spikes)
  smoothed = cell(size(lfp_pow_spikes{i}));
  for j=1:numel(lfp_pow_spikes{i})
    tmp = cellfun(@(c) calculateSmoothedSpikes(c,gausswin(70),select_range(1),diff(select_range)),...
                 lfp_pow_spikes{i}{j}, 'UniformOutput', false);
    smoothed{i}{j} = cellArray2mat(tmp);
  end
  plt = plotMeanAndStds(cellArray2mat(smoothed{i}),'x', select_range(1):select_range(2)-1);
  saveFigures(plt, fullfile(plot_save_dir, 'all_days/spikes', sprintf('psth_shape_%s.png', lfp_pow_names{i})));
end
%% ISI
isi = cell(1,numel(dhs));
for i=1:numel(spikes)
  isi{i} = cell(1,numel(spikes{i}));
  for j=1:numel(spikes{i})
    isi{i}{j} = cellArray2mat(cellfun(@diff, spikes{i}{j}, 'UniformOutput', false));
  end
  isi{i} = cellArray2mat(isi{i});
end
isi = cellArray2mat(isi);
plt = figure('Visible', 'off');
histogram(isi);
saveFigures(plt, fullfile(plot_save_dir, 'all_days/spikes', 'isi_all.png'));
%% ISI by success/failure
tr_isi = cell(size(tr_spikes));
for i=1:numel(tr_spikes)
  tr_isi{i} = cell(size(tr_spikes{i}));
  for j=1:numel(tr_spikes{i})
    for k=1:numel(tr_spikes{i}{j})
      tr_isi{i}{j}{k} = cellArray2mat(cellfun(@diff, tr_spikes{i}{j}{k}, 'UniformOutput', false));
    end
    tr_isi{i}{j} = cellArray2mat(tr_isi{i}{j});
  end
  tr_isi{i} = cellArray2mat(tr_isi{i});
  plt = figure('Visible', 'off');
  histogram(tr_isi{i})
  saveFigures(plt, fullfile(plot_save_dir, 'all_days/spikes', sprintf('isi_%s.png', names{i})));
end
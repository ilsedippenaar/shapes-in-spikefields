function [split_spikes, split_lfps] = splitByPowLabels(spikes, lfps, trial_pow_labs, num_splits, freq_idx)
split_spikes = cell(num_splits,2); % low/high power x S/F
split_lfps = cell(size(split_spikes));
for i=[1,2] % S/F
  for j=1:num_splits % LFP power levels (i.e. low/medium/high)
    split_spikes{j,i} = cell(size(spikes{i}));
    split_lfps{j,i} = cell(size(spikes{i}));
    for k=1:numel(spikes{i}) % day
      split_spikes{j,i}{k} = splitUnitsOrElectrodes(spikes{i}{k}, i, j, k);
      split_lfps{j,i}{k} = splitUnitsOrElectrodes(lfps{i}{k}, i, j, k);
    end
  end
end
  function split = splitUnitsOrElectrodes(data, s_f, level, day)
    % params: success/failure, low/high, day
    split = cell(size(data)); 
    for l=1:numel(data) % number units or electrodes
      if ~isempty(data{l})
        % spikes have 1 x trials cell array 
        % lfps have time x trials matrix
        level_labels = trial_pow_labs{s_f}{day}(freq_idx,:);
        split{l} = data{l}(:,level_labels == level);
      end
    end
  end
end
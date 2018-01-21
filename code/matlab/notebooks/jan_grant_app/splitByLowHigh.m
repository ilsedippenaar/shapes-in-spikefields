function [split_spikes, split_lfps] = splitByLowHigh(spikes, lfps, is_high_tr, freq_idx)
split_spikes = cell(2,2); % low/high power x S/F
split_lfps = cell(size(split_spikes));
for i=[1,2] % S/F
  for j=[1,2] % LFP power (low/high)
    split_spikes{j,i} = cell(size(spikes{i}));
    split_lfps{j,i} = cell(size(spikes{i}));
    for k=1:numel(spikes{i}) % day
      split_spikes{j,i}{k} = splitUnitsOrElectrodes(spikes{i}{k}, i, j, k);
      split_lfps{j,i}{k} = splitUnitsOrElectrodes(lfps{i}{k}, i, j, k);
    end
  end
end
  function split = splitUnitsOrElectrodes(data, s_f, l_h, day)
    % params: success/failure, low/high, day
    split = cell(size(data)); 
    for l=1:numel(data) % number units or electrodes
      if ~isempty(data{l})
        % spikes have 1 x trials cell array 
        % lfps have time x trials matrix
        high_idxs = is_high_tr{s_f}{day}(freq_idx,:);
        if l_h == 1 % l_h = 1 --> low power, so negate is_high_tr
          high_idxs = ~high_idxs;
        end
        split{l} = data{l}(:,high_idxs);
      end
    end
  end
end
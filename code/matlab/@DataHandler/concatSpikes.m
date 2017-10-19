function continuous_spikes = concatSpikes(obj, trial_struct, valid_unit_selec)
% TODO: use cellArray2mat here
% calculate sizes for preallocation
spike_lens = zeros(obj.num_trials, obj.num_units);
for i = 1:obj.num_trials
  tr = trial_struct(i);
  spike_lens(i, :) = cellfun(@numel, tr.sortsp(valid_unit_selec)) ...
    + cellfun(@numel, tr.intsort(valid_unit_selec));
end
cumm_spike_lens = cumsum(spike_lens, 1); % sum down each column
continuous_spikes = cell(1, obj.num_units);
for i = 1:obj.num_units
  continuous_spikes{i} = zeros(1, cumm_spike_lens(end, i));
end

% now just set the spikes according to the indices calculated above
cumm_spike_lens = [zeros(1, obj.num_units); cumm_spike_lens];
unit_nums = find(valid_unit_selec);
for i=1:obj.num_units
  for j=1:obj.num_trials
    tr = trial_struct(j);
    idxs_to_set = (cumm_spike_lens(j,i) + 1) : cumm_spike_lens(j+1,i);
    continuous_spikes{i}(idxs_to_set) = [tr.sortsp{unit_nums(i)}; ...
                                         tr.intsort{unit_nums(i)}];
  end
end

% sort data and exclude duplicates
for i=1:obj.num_units
  continuous_spikes{i} = unique(continuous_spikes{i});
  if ~issorted(continuous_spikes{i})
    warning('Spikes not sorted');
    continuous_spikes{i} = sort(continuous_spikes{i});
  end
end
end
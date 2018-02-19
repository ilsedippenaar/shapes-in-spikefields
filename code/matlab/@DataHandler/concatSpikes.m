function cont_spikes = concatSpikes(obj, trial_struct, valid_unit_selec)
if ~all(isfield(trial_struct, {'sortsp', 'intsort'}))
  cont_spikes = [];
  return
end
% calculate sizes for preallocation
cont_spikes = cell(1, obj.num_trials * 2);
cont_spikes(1:2:end) = {trial_struct.sortsp};
cont_spikes(2:2:end) = {trial_struct.intsort};
for i=1:numel(cont_spikes)
  cont_spikes{i} = cont_spikes{i}(valid_unit_selec);
end
cont_spikes = combineCellArrays('single', cont_spikes{:});
cont_spikes = cellfun(@(c) c(~isnan(c)), cont_spikes, 'UniformOutput', false);

% sort data and exclude duplicates
for i=1:numel(cont_spikes)
  cont_spikes{i} = unique(cont_spikes{i});
  if ~issorted(cont_spikes{i})
    warning('Spikes not sorted');
    cont_spikes{i} = sort(cont_spikes{i});
  end
end
end
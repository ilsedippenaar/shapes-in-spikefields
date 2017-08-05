function dhs = split(obj, betweens)
assert(isnumeric(betweens) && size(betweens,1)==2);

% Initialize
trial_times = [obj.trials.sections];
trial_starts = [trial_times{1:obj.num_trial_sections+1:end}];
trial_ends = [trial_times{obj.num_trial_sections+1:obj.num_trial_sections+1:end}];
n = size(betweens,2);
names = fieldnames(obj);
dhs = cell(1,n);

for i=1:n
  dh = [];
  for j=1:numel(names)
    dh.(names{j}) = obj.(names{j});
  end
  
  % To be changed:
  % fixate, noise, shape, saccade
  % trials, spikes, lfps
  % num_trials
  % number_on_date
  for name={'fixate', 'noise', 'shape', 'saccade'}
    times = dh.(name{1});
    dh.(name{1}) = times(and(times >= betweens(1,i), times < betweens(2,i))) - betweens(1,i) + 1;
  end
  
  start_trial_num = binarySearch(trial_starts, betweens(1,i), ']', true);
  end_trial_num = binarySearch(trial_ends, betweens(2,i), '(', true);
  dh.trials = obj.trials(start_trial_num:end_trial_num);
  dh.num_trials = numel(dh.trials);
  for j=1:dh.num_trials
    dh.trials(j).saccade = dh.trials(j).saccade - betweens(1,i) + 1;
    for k=1:dh.num_trial_sections+1
      dh.trials(j).sections{k} = dh.trials(j).sections{k} - betweens(1,i) + 1;
    end
  end
  
  dh.lfps = cellArray2mat(obj.select('type', 'lfp', 'between', betweens(:,i), 'melt', true), 'int16');
  dh.spikes = obj.select('type', 'spike', 'between', betweens(:,i), 'melt', true);
  for j=1:obj.num_units
    dh.spikes{j} = dh.spikes{j} - betweens(1,i) + 1;
  end
  
  dh.number_on_date = i;
  dhs{i} = DataHandler('dh_struct', dh, 'clean', false);
end
end
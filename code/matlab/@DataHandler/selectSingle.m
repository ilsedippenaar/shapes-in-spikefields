function out = selectSingle(obj, type, between, buffer, index, trial_num, trial_section)
out = [];
switch trial_section
  case -1
    start = obj.trials(trial_num).sections{5};
    if isempty(start)
      error('There must a a shape for selecting between shape and saccade');
    end
    section = [start obj.trials(trial_num).saccade];
    if numel(section) ~= 2
      error('There must a a saccade for selecting between shape and saccade');
    end
  case -2
    start = obj.trials(trial_num).saccade;
    if isempty(start)
      error('There must a a saccade for selecting between saccade and stop');
    end
    section = [start obj.trials(trial_num).sections{6}];
  otherwise
    start = obj.trials(trial_num).sections{trial_section};
    if isempty(start)
      return
    end
    section = [start min([obj.trials(trial_num).sections{trial_section+1:end}])]; % this should always be in range
end
assert(section(2) >= section(1));
section = section + [-buffer(1) buffer(2)];
times = intersectInterval(between, section);
if isempty(times)
  return
end

if strcmpi(type, 'spike')
  spikes = obj.spikes{index};
  start = binarySearch(spikes, times(1), false);
  stop  = binarySearch(spikes, times(2), true);
  if ~isempty(start) && ~isempty(stop)
    out = spikes(start:stop);
  end
elseif strcmpi(type, 'lfp')
  start = max(1, floor(times(1)));
  stop = min(size(obj.lfps, 1), ceil(times(2))-1);
  out = obj.lfps(start:stop, index);
end
end
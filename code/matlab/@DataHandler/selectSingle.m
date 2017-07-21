function out = selectSingle(obj, type, between, buffer, index, trial_num, trial_section)
out = [];
start = obj.trials(trial_num).sections{trial_section};
if isempty(start)
  return
end
section = [start min([obj.trials(trial_num).sections{trial_section+1:end}])]; % this should always be in range
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
  stop = min(size(obj.lfps, 1), ceil(times(2)));
  out = obj.lfps(start:stop, index);
end
end
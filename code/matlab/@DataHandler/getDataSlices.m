function [out,times] = getDataSlices(obj, type, select_range, conditions, varargin)
% assumes conditons.vec covers of the entire data time period (i.e. that
% all data recorded are reflective of the conditions specified and that
% there aren't data recorded when conditions were not).
p = inputParser;
p.addParameter('data_type', []);
p.parse(varargin{:});
args = p.Results;

% validate and parse input
if strcmp(type, 'lfp')
  data = obj.lfps;
  data_type = class(obj.lfps);
else
  data = obj.spikes;
  data_type = class(obj.spikes{1});
end
if ~isempty(args.data_type)
  data_type = args.data_type;
end
start_time = 0;
if strcmp(type, 'lfp')
  end_time = size(data,1) + start_time - 1;
else
  end_time = min(cellfun(@(c) c(end), data(~cellfun(@isempty, data))));
end

% make difference vectors for the negated conditions
for i=1:numel(conditions)
  if ischar(conditions(i).vec)
    conditions(i).vec = obj.(conditions(i).vec);
  end
  if conditions(i).negate
    % this is a bit of a hack to ensure none conditions lie within the
    % allowed range
    idx = binarySearch(conditions(i).vec, start_time, '(');
    if isempty(idx)
      idx = 0;
    end
    conditions(i).vec = [conditions(i).vec(1:idx) start_time conditions(i).vec(idx+1:end)];
    
    idx = binarySearch(conditions(i).vec, end_time, '(');
    if isempty(idx)
      idx = numel(conditions(i).vec);
    end
    conditions(i).vec = [conditions(i).vec(1:idx) end_time conditions(i).vec(idx+1:end)];
    
    conditions(i).diff = diff(conditions(i).vec);
  end
end

% initialize
time = max(start_time, start_time - select_range(1)); % have select_range(1) data points before the "0" time

% get all valid times
times = [];
while time + select_range(2) - 1 <= end_time
  i = 0;
  j = 1;
  while i ~= numel(conditions)
    next_time = getNextValidTime(conditions(j), time);
    if isempty(next_time)
      time = [];
      break
    end
    if next_time == time
      i = i + 1;
    else
      time = next_time;
      i = 1;
    end
    j = mod(j, numel(conditions)) + 1;
  end
  if isempty(time) || time + select_range(2) - 1 > end_time
    break
  end
  times(end+1) = time;
  time = time + diff(select_range); % simulate copy
end

% copy data
if strcmp(type, 'lfp')
  out = cell(1,size(data,2));
else
  out = cell(1,numel(data));
end
convertFunc = str2func(data_type);
if strcmp(type, 'lfp')
  idxs = times - start_time + 1;
  for i=1:size(data,2)
    out{i} = zeros(diff(select_range), numel(idxs), data_type);
    for j=1:numel(idxs)
      out{i}(:,j) = convertFunc(data(idxs(j)+select_range(1) : idxs(j)+select_range(2)-1,i));
    end
  end
else
  for i=1:numel(data)
    out{i} = cell(1,numel(times));
    for j=1:numel(times)
      start_idx = binarySearch(data{i}, times(j)+select_range(1), ']');
      stop_idx = binarySearch(data{i}, times(j)+select_range(2), '(');
      out{i}{j} = convertFunc(data{i}(start_idx:stop_idx) - double(times(j)));
    end
  end
end

  function next_time = getNextValidTime(condition, current_time)
    % can (and should sometimes) return current_time
    if condition.negate
      if hasElementBetween(condition.vec, current_time+condition.range(1), current_time+condition.range(2))
        diff_idx = binarySearch(condition.vec, current_time + condition.range(1), ']', true);
        next_diff_idx = [];
        for i_=diff_idx:numel(condition.diff)
          if condition.diff(i_) > diff(condition.range)
            next_diff_idx = i_;
            break
          end
        end
        if isempty(next_diff_idx)
          next_time = []; % off the end of the array
        else
          next_time = condition.vec(next_diff_idx)+1-condition.range(1); % finds the earliest time with large enough diff
        end
      else
        next_time = current_time; % the current time works
      end
    else
      next_idx = binarySearch(condition.vec, current_time + condition.range(1), ']', true);
      if isempty(next_idx)
        next_time = [];
      else
        next_time = max(current_time, condition.vec(next_idx) - condition.range(2) + 1);
      end
    end
  end
end
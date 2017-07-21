function out = getLfpSlices(obj, spec, start_time, end_time, varargin)
p = inputParser;
p.addParameter('selection_params', {}, @iscell);
p.addParameter('lfps', [], @(x) isnumeric(x) || (iscell(x) && all(cellfun(@isnumeric, x))));
p.addParameter('data_type', 'int16');
p.parse(varargin{:});
args = p.Results;

if isempty(args.selection_params) && isempty(args.lfps)
  lfps = obj.select('type', 'lfp', 'between', [start_time end_time], 'melt', true);
else
  if ~isempty(args.selection_params)
    lfps = obj.select(args.selection_params{:});
  elseif isnumeric(args.lfps)
    lfps = cell(1,size(args.lfps,2));
    for i=1:numel(lfps)
      lfps{i} = args.lfps(:,i);
    end
  else
    lfps = args.lfps; % this should be a cell array
  end
end

n = numel(lfps{1});
all_events = [obj.trials.sections];
fixate = [all_events{3:7:end}] - start_time + 1;
fixate = int32(fixate(fixate>0));
noise = [all_events{4:7:end}] - start_time + 1;
noise = int32(noise(noise>0));
shape = [all_events{5:7:end}] - start_time + 1;
shape = int32(shape(shape>0));
saccade = [obj.trials.saccade] - start_time + 1;
saccade = int32(saccade(saccade>0));
all_events = sort([fixate noise shape saccade]); % for the none condition

% parse spec
params = strsplit(spec, ',', 'CollapseDelimiters', false);
buffer = zeros(1,2, 'int32');
if ~isempty(params{1})
  buffer(1) = int32(str2double(params{1}));
end
if ~isempty(params{end})
  buffer(2) = int32(str2double(params{end}));
end

cond_names = [string(params{2}) string(params{5})];
cond_sections = zeros(1,2, 'int32');
cond_vecs = cell(1,2);
for i=1:2
  switch cond_names(i)
    case 'noise'
      cond_vecs{i} = noise;
    case 'shape'
      cond_vecs{i} = shape;
    case 'fixate'
      cond_vecs{i} = fixate;
    case 'saccade'
      cond_vecs{i} = saccade;
    case 'none'
      cond_vecs{i} = [];
    case 'any'
      cond_vecs{i} = [];
    otherwise
      error('Unsupported condition name: %s', cond_names(i));
  end
  if cond_names(i) ~= 'any'
    if i == 1
      num = 3;
    else
      num = 6;
    end
    cond_sections(i) = int32(str2double(params{num}));
  end
end

section = 0;
if ~isempty(params{4})
  section = int32(str2double(params{4}));
end

% this is an optimization to speed up what would otherwise be very slow
sizes = cellfun(@numel, cond_vecs) .* double(cond_sections);
tmp_min = min(sizes(and(sizes > 0, cond_names ~= 'any')));
if isempty(tmp_min)
  align_num = [];
else
  align_num = find(sizes == tmp_min);
end

% allocate space
instances = 0;
runFuncInValidSections(@countNumValid);
out = cell(1,numel(lfps));
for idx=1:numel(lfps)
  out{idx} = zeros(section + sum(buffer), instances, args.data_type);
end

% copy lfps
curr_instance = 1;
convert_func = str2func(args.data_type);
runFuncInValidSections(@copyLfps);

  function countNumValid(varargin)
    instances = instances + 1;
  end
  function copyLfps(idx)
    for j=1:numel(lfps)
      out{j}(:,curr_instance) = convert_func(lfps{j}(idx-buffer(1):idx+buffer(2)+section-1));
    end
    curr_instance = curr_instance + 1;
  end


  function runFuncInValidSections(f)
    % setting some parameters
    i = max(buffer(1), cond_sections(1)) + 1;
    if ~isempty(align_num)
      if align_num == 1
        i_offset = 1;
      else
        i_offset = -cond_sections(align_num) - section + 1;
      end
      num_in_section = cond_sections(align_num);
    else
      num_in_section = 1;
    end
    
    % initializing i for aligning
    if ~isempty(align_num)
      if align_num == 1
        lower_bound = i - cond_sections(align_num);
      else
        lower_bound = i + section;
      end
      next_valid = find(cond_vecs{align_num} >= lower_bound, 1);
      if isempty(next_valid)
        return
      end
      i = cond_vecs{align_num}(next_valid) + i_offset;
    end
    
    while i <= n - section - max(buffer(2), cond_sections(2))
      for j=1:num_in_section
        [valid1, none1] = checkValid(cond_names(1), cond_vecs{1}, i-cond_sections(1), i);
        [valid2, none2] = checkValid(cond_names(2), cond_vecs{2}, i+section, i+section+cond_sections(2));
        if valid1 && valid2
          f(i);
          i = i + section + sum(buffer);
          break % out of for
        elseif ~isempty(none1) || ~isempty(none2)
          if ~isempty(none1)
            i_for_none1 = none1 + cond_sections(1);
          else
            i_for_none1 = i;
          end
          if ~isempty(none2)
            i_for_none2 = none2 - section;
          else
            i_for_none2 = i;
          end
          i = max(i_for_none1, i_for_none2);
          break % out of for
        else
          i = i + 1;
        end
      end
      if ~isempty(align_num)
        if align_num == 1
          lower_bound = i - cond_sections(align_num);
        else
          lower_bound = i + section;
        end
        next_valid = find(cond_vecs{align_num} >= lower_bound, 1);
        if isempty(next_valid)
          return
        end
        i = cond_vecs{align_num}(next_valid) + i_offset;
      end
    end
  end

  function [valid, none_time] = checkValid(condition, cond_vec, start, stop)
    none_time = [];
    if condition == 'none'
      valid = ~hasElementBetween(all_events, start, stop);
      if ~valid
        none_time = all_events(binarySearch(all_events, stop, '(')) + 1;
      end
    elseif condition == 'any'
      valid = true;
    else
      valid = hasElementBetween(cond_vec, start, stop);
    end
  end
end
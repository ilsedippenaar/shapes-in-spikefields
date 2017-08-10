function out = select(obj, varargin)
%SELECT Select portions of the data matching a few parameters.
%
%   type: char, default: 'spike', one of {'spike', 'lfp'}
%   melt: logical, default: false
%     If true, trial_result and trial_section are ignored and the data are
%     concatenated together into a single vector per electrode
%   units: char, default: 'ms', one of {'ms', timestamp', 'seconds'}
%     Indicates the units for the between and buffer parameters
%   between: 1x2 double array, default: [0, Inf]
%     Selects data only between the time points given.
%   buffer: 1x2 double array, default: [0, 0]
%     Includes data before and after the selection (if possible)
%   electrode_nums, 1xn double array, default: 1:DataHandler.num_lfp_electrodes
%   unit_nums, 1xn double array, default: 1:DataHandler.num_units
%
%   trial_nums: 1xn double array, default: 1:numel(trials)
%   trial_result: 1xn cell array, default: {'true_positive', 'true_negative', 'false_positive', 'false_negative' 'failed'}
%     The cell array can be any size and can define trial results in
%     combination (i.e. include both true_positive and true_negative
%     outcomes).
%   trial_section: 1xn double array or char or cell array of char, default: 2:6
%     If trial_section is a char or cell array of char, then it must be one
%     of the keys in DataHandler.section_map. The section coding scheme for
%     numbered trial sections follows the diagram below.
%
%        1  | 2  |  3 |   4   |   5   |   6
%           |    |    |       |       |     
%         start fix  noise  shape    end
%
%     For example, if trial_section == 2, then only the sections of trials
%     between the start of the trial and a fix event would be returned. If 
%     an event was not defined for that trial (i.e. no saccade occurred), 
%     then data in the time range which would have included it is returned. 
%     So, if trial_result == {'true_negative'} and trial_section == 4, then 
%     data between the fix event and the end of the trial is returned. A 
%     value of 1 indicates the time between the end of the previous trial 
%     and the beginning on the current one, and 6 indicates the time 
%     between the end of this trial and the beginning on the next one.

p = inputParser;
stringIn = @(string, set) any(strcmpi(string, set));
default_trial_result = {'true_positive', 'true_negative', 'false_positive', 'false_negative' 'failed'};

p.addParameter('type', 'spike', @(x) stringIn(x, {'spike', 'lfp'}));
p.addParameter('melt', false, @islogical);
p.addParameter('units', 'ms', @(x) stringIn(x, {'ms', 'timestamp', 'seconds'}));
p.addParameter('between', [0 Inf], @(x) numel(x) == 2);
p.addParameter('electrode_nums', 1:obj.num_lfp_electrodes);
p.addParameter('unit_nums', 1:obj.num_units);
p.addParameter('buffer', [0 0], @(x) numel(x) == 2);
p.addParameter('trial_nums', 1:obj.num_trials);
p.addParameter('trial_result', default_trial_result, @(x) all(cellfun(@(s) stringIn(s, default_trial_result), x)));
p.addParameter('trial_section', 2:6);

p.parse(varargin{:});
args = p.Results;

% convert buffer and between to ms
switch args.units
  case 'timestamp'
    between = args.between ./ obj.ts_per_ms;
    buffer = args.buffer ./ obj.ts_per_ms;
  case 'seconds'
    between = args.between .* 1e3;
    buffer = args.buffer .* 1e3;
  otherwise
    between = args.between;
    buffer = args.buffer;
end
between = [floor(between(1)) ceil(between(2))];
buffer = [floor(buffer(1)) ceil(buffer(2))];

if iscell(args.trial_section)
  trial_sections = cellfun(@obj.section_map, args.trial_section);
elseif ischar(args.trial_section)
  trial_sections = obj.section_map(args.trial_section);
else
  trial_sections = args.trial_section;
  assert(min(trial_sections) >= -2 && min(trial_sections) ~= 0);
  assert(max(trial_sections) <= obj.num_trial_sections);
end
      
% get possible trial nums from trial_nums and trial_result
if args.melt
  trial_nums = args.trial_nums;
  start_trial = floor(obj.trials(trial_nums(1)).sections{1}); % assume trial_nums is contiguous
  stop_trial = ceil(obj.trials(trial_nums(end)).sections{end});
  between = intersectInterval(between, [start_trial stop_trial]);
else
  all_results = {obj.trials.result};
  trial_nums = find(cellfun(@(c) stringIn(c, args.trial_result), all_results));
  trial_nums = intersect(args.trial_nums, trial_nums);
end

switch lower(args.type)
  case 'spike'
    indices = args.unit_nums;
  case 'lfp'
    indices = args.electrode_nums;
end

out = cell(1, numel(indices));
for i=1:numel(indices)
  index = indices(i);
  if args.melt
    switch lower(args.type)
      case 'spike'
        spikes = obj.spikes{index};
        % only get inclusive values in range
        start = binarySearch(spikes, between(1), ']', true);
        stop = binarySearch(spikes, between(2), '(', true);
        out{i} = spikes(start:stop);
      case 'lfp'
        lfp = obj.lfps(:,index);
        start = max(1, between(1));
        stop = min(numel(lfp), between(2));
        out{i} = lfp(floor(start):ceil(stop)-1);
    end
  else
    out{i} = cell(numel(trial_nums), numel(trial_sections));
    for j=1:numel(trial_nums)
      trial_num = trial_nums(j);
      for k=1:numel(trial_sections)
        trial_section = trial_sections(k);
        out{i}{j, k} = obj.selectSingle(args.type, between, buffer, ...
          index, trial_num, trial_section);
      end
    end
  end
end
end
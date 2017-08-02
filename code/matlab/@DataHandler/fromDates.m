function data_handlers = fromDates(dates, trial_dir, lfp_dir, data_handler_dir, varargin)
if ischar(dates)
  num = 1;
else
  num = numel(dates);
end

valid(num).files = [];
valid(num).date = [];
for i=1:num
  if iscell(dates)
    d = dates{i};
  elseif ischar(dates)
    d = dates;
  else % datetime objects
    d = dates(i);
  end
  date_string = datestr(d, 'yyyy-mm-dd');
  dh_file = fullfile(data_handler_dir, [date_string '.mat']);
  if exist(dh_file, 'file') == 2
    valid(i).files = dh_file;
    valid(i).date = date_string;
  else
    trial_files = dir(fullfile(trial_dir, [date_string '*']));
    lfp_files = dir(fullfile(lfp_dir, [date_string '*']));
    prefixes = {'More than one', 'No'};
    data_names = {'trial', 'lfp'};
    files = {{trial_files.name} {lfp_files.name}};
    is_valid = true;
    for j=1:2
      if numel(files{j}) > 1
        fprintf('%s %s file found for %s\n', prefixes{1}, data_names{j}, date_string);
        is_valid = false;
      elseif numel(files{j}) == 0
        fprintf('%s %s file found for %s\n', prefixes{2}, data_names{j}, date_string);
        is_valid = false;
      end
    end
    if is_valid
      valid(i).files = files;
      valid(i).date = date_string;
    else
      fprintf('\n');
    end
  end
end

valid_idxs = find(~cellfun(@isempty, {valid.files}));
total_valid = numel(valid_idxs);
if num > 1 % i.e. more than one date requested
  data_handlers(total_valid).dh = [];
end

if total_valid == 0
  data_handlers = [];
  return
end

fprintf('Found %d valid dates, loading now...\n', total_valid);
for i=1:total_valid
  valid_idx = valid_idxs(i);
  dh_file = fullfile(data_handler_dir, [valid(valid_idx).date '.mat']);
  if numel(valid(valid_idx).files) > 1 && iscell(valid(valid_idx).files)
    fprintf('Making new DataHandler for %s\n', valid(valid_idx).date);
    % overwrite save name if it exists
    p = inputParser;
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    args = p.Unmatched;
    args.save_name = dh_file;
    args.date = valid(valid_idx).date;
    dh = DataHandler.fromFile(fullfile(trial_dir, valid(valid_idx).files{1}{1}), ...
                              fullfile(lfp_dir, valid(valid_idx).files{2}{1}), args);
  else
    fprintf('Loading DataHandler from cached state for %s\n', valid(valid_idx).date);
    tmp = load(dh_file);
    dh = tmp.dh;
  end
  
  if num > 1
    data_handlers(i).dh = dh;
  else
    data_handlers = dh; % at most 1 requested
  end
end
if num > 1
  data_handlers = [data_handlers.dh];
end
end
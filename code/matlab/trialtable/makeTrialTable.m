function tbl = makeTrialTable(days, trial_dir, lfp_dir, dh_dir, save_dir, select_type, varargin)
if ischar(days)
  days = {days};
end
tbls = cell(1,numel(days));
for i=1:numel(days)
  if iscell(days)
    d = days{i};
  else % datetime objects
    d = days(i);
  end
  day_str = datestr(d, 'yyyy-mm-dd');
  save_name = fullfile(save_dir, day_str);
  if exist([save_name '.mat'], 'file') == 2
    if nargout > 0
      tmp = load(save_name);
      tbls{i} = tmp.tbl;
    end
  else
    dh = DataHandler.fromDates(d, trial_dir, lfp_dir, dh_dir, varargin{:});
    if isempty(dh)
      continue
    end
    dhs = dh.split();
    dhs = [dhs{:}];
    dhs = dhs([dhs.num_trials] ~= 0);
    dhs = dhs(~cellfun(@isempty, {dhs.lfps}));
    tbls{i} = arrayfun(@(x) x.toTrialTable, dhs, 'UniformOutput', false);
    tbls{i} = vertcat(tbls{i}{:});
    if ~isempty(select_type)
      tbls{i} = select(tbls{i}, select_type);
    end
    if ~isempty(save_dir)
      tmp = [];
      tmp.tbl = tbls{i};
      save(save_name, '-struct', 'tmp');
    end
  end
end
tbl = vertcat(tbls{:});
end
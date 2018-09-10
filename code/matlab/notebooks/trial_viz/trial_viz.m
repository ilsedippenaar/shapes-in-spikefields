monkey_name = 'zorin';
folder_name = 'trial_viz';

trial_dir = fullfile(data_dir, 'trials', monkey_name);
lfp_dir = fullfile(data_dir, 'lfps', monkey_name);
dh_dir = fullfile(cache_dir, 'data_handlers', monkey_name);
save_dir = fullfile(cache_dir, 'trial_tables', 'all', monkey_name);
python_save_dir = fullfile(data_dir, '../python', folder_name);

if strcmp(monkey_name, 'jaws')
  dates = datetime('2013-05-1'):datetime('2013-07-19');
  learn_start = datetime('2013-06-3');
else
  dates = datetime('2011-7-20'):datetime('2011-12-20');
  learn_start = datetime('2011-11-7');
end
%% 
tbl = makeTrialTable(dates, trial_dir, lfp_dir, dh_dir, [], [], ...
  'split_dh', false, ...
  'dh_args', {...
    'clean', true, ...
    'min_reaction_time', params.length_min_reaction_time, ...
    'monkey_name', monkey_name});
%tbl = tbl(~strcmp(tbl.result, 'failed'),:);
%tbl = tbl(~isnan(tbl.shape),:);
tbl.trace = cellfun(@(c) single(mean(c,2)), tbl.lfps, 'UniformOutput', false);
grouped = groupBy(tbl, 'date');
%%
for i=1:numel(grouped)
  ss = [];
  ss.num = int32(grouped{i}.num);
  ss.lfp = grouped{i}.trace;
  spikes = grouped{i}.spikes;
  for j=1:numel(spikes)
    spikes{j} = cellfun(@(c) (c-grouped{i}.start(j))/1000, spikes{j}, 'UniformOutput', false);
  end
  ss.spikes = spikes;
  ss.events = struct('fixate', single(grouped{i}.fixate - grouped{i}.start)/1000, ...
                      'noise', single(grouped{i}.noise - grouped{i}.start)/1000, ...
                      'shape', single(grouped{i}.shape - grouped{i}.start)/1000, ...
                      'saccade', single(grouped{i}.saccade - grouped{i}.start)/1000);
  ss.shapecoh = single(grouped{i}.shapecoh);
  ss.result = grouped{i}.result;
  date = grouped{i}.date{1};
  name = grouped{i}.monkey_name{1};
  save(fullfile(python_save_dir, sprintf('%s_%s', date, name)), '-struct', 'ss');
end
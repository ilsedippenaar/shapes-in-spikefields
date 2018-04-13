function tbl = toTrialTable(obj)
% columns = id, date, num, monkey_name, lfps, spikes, electrode_mapping
%           start, fixate, noise, shape, saccade, stop, result

% sometimes there are 2 time points in a condition - these are just useless
% trials
valid_trials = arrayfun(@(x) all(cellfun(@numel, x.sections) < 2), obj.trials);
trials = obj.trials(valid_trials);
sections = vertcat(trials.sections);
saccade = {trials.saccade}';
sections = [sections(:,2:5), saccade, sections(:,6)]; % start, fix, noise, shape, saccade, stop
for i=1:numel(sections)
  if isempty(sections{i})
    sections{i} = nan;
  end
end
disp(obj.date)
sections = cellfun(@(c) c, sections); % num_trials x 6 double
% sanity check - remove rows that have nan for either start or stop
sections = sections(~or(isnan(sections(:,1)), isnan(sections(:,end))),:);

num_trials = size(sections,1);
[ids, lfps, spikes] = deal(cell(num_trials, 1));
mat = @(x) [x{:}];
for i=1:num_trials
  ids{i} = [obj.date, '_', num2str(obj.number_on_date), '_', num2str(i)];
  lfps{i} = mat(obj.select('type', 'lfp', 'between', sections(i, [1 end]), 'melt', true));
  spikes{i} = obj.select('type', 'spike', 'between', sections(i, [1 end]), 'melt', true);
end

sections = mat2cell(sections, size(sections,1), ones(1,size(sections,2)));
date = repelem({obj.date}, num_trials, 1);
num = (1:num_trials)';
monkey_name = repelem({obj.monkey_name}, num_trials, 1);
electrode_mapping = repelem({obj.electrode_mapping},num_trials,1);
shapeid = [trials.shapeid]';
shapecoh = [trials.shapecoh]';
result = {trials.result}';
tbl = table(ids, date, num, monkey_name, lfps, spikes, electrode_mapping, sections{:}, shapeid, shapecoh, result, ...
  'VariableNames', {'id','date','num','monkey_name', 'lfps', 'spikes','electrode_mapping', ...
                    'start', 'fixate', 'noise', 'shape', 'saccade', 'stop', 'shapeid', 'shapecoh', 'result'});
end
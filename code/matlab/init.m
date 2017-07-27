%%
if exist('dh', 'var') ~= 1
  dh = DataHandler.fromDates('2013-05-1', ...
    fullfile(data_dir, 'trials'), ...
    fullfile(data_dir, 'lfps'), ...
    fullfile(data_dir, 'data_handlers'));
end
%% Subsetting data (for 05-01)
fprintf('Selecting data\n');
start_time = (24*60 + 10)*1000;
lfps = dh.select('between', [start_time, size(dh.lfps,1)-30*1000], 'melt', true, 'type', 'lfp');
spikes = dh.select('between', [start_time, size(dh.lfps,1)-30*1000], 'melt', true, 'type', 'spike');
%% Making LFP slices from conditions
fprintf('Making slices\n');
select_range = [0 400];
fixate_cond = cell2struct([...
  {'fixate', false, [0 1]}; ...
  {'shape', true, select_range}; ...
  {'noise', true, select_range}; ...
  {'saccade', true, select_range}], {'vec', 'negate', 'range'}, 2);
shape_cond = cell2struct([...
  {'fixate', true, select_range}; ...
  {'shape', false, [0 1]}; ...
  {'noise', true, select_range}; ...
  {'saccade', true, select_range}], {'vec', 'negate', 'range'}, 2);
noise_cond = cell2struct([...
  {'fixate', true, select_range}; ...
  {'shape', true, select_range}; ...
  {'noise', false, [0 1]}; ...
  {'saccade', true, select_range}], {'vec', 'negate', 'range'}, 2);
saccade_cond = cell2struct([...
  {'fixate', true, select_range}; ...
  {'shape', true, select_range}; ...
  {'noise', true, select_range}; ...
  {'saccade', false, [0 1]}], {'vec', 'negate', 'range'}, 2);
null_cond = cell2struct([...
  {'fixate', true, select_range}; ...
  {'shape', true, select_range}; ...
  {'noise', true, select_range}; ...
  {'saccade', true, select_range}], {'vec', 'negate', 'range'}, 2);
conditions = {fixate_cond, shape_cond, noise_cond, saccade_cond, null_cond};
names = {'fixate', 'shape', 'noise', 'saccade', 'null'};
lfp_slices = cell(1,numel(names));
for i=1:numel(lfp_slices)
  lfp_slices{i} = dh.getDataSlices(lfps, 'lfp', select_range, conditions{i}, 'start_time', start_time, 'data_type', 'single');
end
clear fixate_cond shape_cond noise_cond saccade_cond null_cond

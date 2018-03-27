params = [];
% Selecting
params.names = {'shape', 'saccade'};
params.length_postnoise_response = 1000;
params.length_presaccade_response = 50;
params.length_min_reaction_time = 150;
% Multitaper
params.T = 0.5;
params.W = 6;
params.window_size = 0.064;
params.Fs = 1000;
% Plotting
params.freq_cutoff = 100;

name = 'jaws';
trial_dir = fullfile(data_dir, 'trials', name);
lfp_dir = fullfile(data_dir, 'lfps', name);
dh_dir = fullfile(cache_dir, 'data_handlers', name);
if strcmp(name, 'jaws')
  %dates = datetime('2013-05-1'):datetime('2013-07-19');
  dates = datetime('2013-05-1'):datetime('2013-05-1');
else
  dates = datetime('2011-07-22'):datetime('2011-12-20');
end

if exist('dhs', 'var') ~= 1
  dhs = DataHandler.fromDates(dates, trial_dir, lfp_dir, dh_dir, ...
    'clean', true, ...
    'min_reaction_time', params.length_min_reaction_time, ...
    'monkey_name', name);
end

all_dhs = cell(1,numel(dhs));
for i=1:numel(dhs)
  all_dhs{i} = dhs(1).split();
  dhs(1) = []; % delete some old dhs to save on memory
end
dhs = [all_dhs{:}];
clear all_dhs
dhs = [dhs{:}];
dhs = dhs([dhs.num_trials] ~= 0);
dhs = dhs(~cellfun(@isempty, {dhs.lfps}));
clear num idx name trial_dir lfp_dir dh_dir dates i j
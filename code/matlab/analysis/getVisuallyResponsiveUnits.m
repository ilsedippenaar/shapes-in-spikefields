function [units,p] = getVisuallyResponsiveUnits(day, monkey_name, min_reaction_time, ...
                                            trial_dir, lfp_dir, dh_dir, all_save_dir)
tbl = makeTrialTable(day, trial_dir, lfp_dir, dh_dir, all_save_dir, [], ...
  'clean', true, ...
  'min_reaction_time', min_reaction_time, ...
  'monkey_name', monkey_name);

win = [-500,500];
valid = [~isnan(tbl.noise), ...
  between(tbl.noise, tbl.fixate - win(1), tbl.saccade - win(2))];
tbl = tbl(all(valid,2),:);
[~,spikes] = getAnalysisData(tbl, win, 'noise');

% just do basic spike counts to see which units give *different* responses
% to a visual stimulus
spikes = vertcat(spikes{:}); % trials x units
if isempty(spikes)
  units = [];
  p = [];
  return
end

no_stim_counts = cellfun(@(c) sum(c < 0), spikes);
stim_counts = cellfun(@numel, spikes) - no_stim_counts;
% one-tailed t-test without assuming variances are equal - gets neurons
% which have greater stim_counts than no_stim_counts (increase firing rate)
[units,p] = ttest2(stim_counts, no_stim_counts, 'Alpha', 0.01, 'Vartype', 'unequal', 'Tail', 'right');
units = logical(units);
end
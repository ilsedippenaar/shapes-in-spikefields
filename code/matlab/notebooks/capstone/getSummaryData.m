function [events, lfp, spikes] = getSummaryData(dh, trial_num)
trial = dh.trials(trial_num);
times = [trial.sections{1}+3000, trial.sections{end}-4000];
cond = [];
cond.vec = times(1);
cond.negate = false;
cond.range = [0,1];

events = [trial.sections{3:5}, trial.saccade]-times(1);

tmp = dh.getDataSlices('lfp', [0, diff(times)], cond);
lfp = mean([tmp{:}],2);

spikes = dh.getDataSlices('spike', [0, diff(times)], cond);
spikes = cellfun(@(c) c{1}, spikes, 'UniformOutput', false);
end
function [con, freqs] = connectivityLFP(lfps, trial_start, align_times)
% lfps = {trial} x time x channels
num_trials = numel(lfps);
trial_len = size(lfps{1},1);
num_chans = size(lfps{1},2);

% 'data' needs fields: trial, sampleinfo, time, fsample, label
data = [];
data.trial = cellfun(@(x) single(x'), lfps, 'UniformOutput', false);
data.sampleinfo = [align_times, align_times+trial_len-1] + trial_start;
data.time = repelem({(trial_start:trial_start+trial_len-1) / 1000}, num_trials, 1);
data.fsample = 1000;
data.label = arrayfun(@(x) ['e' num2str(x)], 1:num_chans, 'UniformOutput', false);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq = 5;
cfg.pad = 'nextpow2';
cfg.output = 'fourier'; % get complex fourier coefficients
cfg.feedback = 'none';

freq = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'ppc';
out = ft_connectivityanalysis(cfg, freq);
con = out.ppcspctrm;
freqs = out.freq;
end
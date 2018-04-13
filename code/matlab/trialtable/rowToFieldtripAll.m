function data = rowToFieldtripAll(lfps, spikes, trial_beg, trial_len, Fs)
if nargin < 5
  Fs = 1000;
end
electrode_labels = strsplit(num2str(1:size(lfps,2)));
electrode_labels = cellfun(@(c) ['e' c], electrode_labels, 'UniformOutput',false);
spike_labels = strsplit(num2str(1:numel(spikes)));
spike_labels = cellfun(@(c) ['u' c], spike_labels, 'UniformOutput', false);

data = [];
data.trial = {[lfps, expandSpikes(spikes, trial_beg, trial_len)]'};
data.time = {(trial_beg:trial_beg+trial_len-1) / Fs};
data.fsample = Fs;
data.label = [electrode_labels, spike_labels]';
end
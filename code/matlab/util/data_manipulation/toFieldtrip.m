function [data, electrode_labels, spike_labels] = toFieldtrip(lfps, spikes, trial_beg, Fs, electrode_label_prefix, spike_label_prefix)
% returns result of ft_appendspike
% lfps = 1 x num_electrodes -> time x num_trials int16/single 
% spikes = 1 x num_units -> 1 x num_trials -> num_spikes x 1 single
% rearrange to have electrodes nested in a trial
if nargin < 6
  spike_label_prefix = 'u';
  if nargin < 5
    electrode_label_prefix = 'e';
  end
end

valid_elecs = ~cellfun(@isempty, lfps);
lfps = lfps(valid_elecs);

if isempty(lfps) || isempty(spikes)
  [data, electrode_labels, spike_labels] = deal([]);
  return
end
num_electrodes = numel(lfps);
num_units = numel(spikes);

lfps = reshape([lfps{:}], [size(lfps{1}), num_electrodes]); % time x trials x electrodes
lfps = permute(lfps, [1,3,2]); % time x electrodes x trials

% 1 x num_trials -> 1 x num_units -> num_spikes x 1 single
spikes = zipCell(spikes{:});

num_trials = size(lfps, 3);
trial_len = size(lfps, 1);
electrode_labels = strsplit(num2str(1:num_electrodes));
electrode_labels = cellfun(@(c) [electrode_label_prefix c], electrode_labels, 'UniformOutput',false);
spike_labels = strsplit(num2str(1:num_units));
spike_labels = cellfun(@(c) [spike_label_prefix c], spike_labels, 'UniformOutput', false);

data = [];
[data.trial,data.time] = deal(cell(1,num_trials));
data.fsample = Fs;
data.label = [electrode_labels, spike_labels]';

for trial=1:num_trials
  expanded_spikes = expandSpikes(spikes{trial}, trial_beg, trial_len);
  data.trial{trial} = double([lfps(:,:,trial), expanded_spikes]');
  data.time{trial} = (trial_beg:trial_beg+trial_len-1) / 1000;
end
end
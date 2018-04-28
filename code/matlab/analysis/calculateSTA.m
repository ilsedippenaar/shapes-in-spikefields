function [sta, time, std_err] = calculateSTA(lfps, spikes, trial_beg, Fs, timwin)
% lfps = days -> electrodes -> time x trials
% spikes = days -> units -> trials -> time
% averaged over all trials, units, and electrodes
cfg = [];
if nargin > 4
  cfg.timwin = timwin; % ft defaults to [-0.1, 0.1]
end
cfg.keeptrials = 'yes';
cfg.latency = 'prestim';

num_days = numel(lfps);
sta = cell(1,num_days);
time = [];
for day=1:num_days
  [ft_data, lfp_labs, spike_labs] = toFieldtrip(lfps{day}, spikes{day}, trial_beg, Fs);
  cfg.channel = lfp_labs;
  num_units = numel(spike_labs);
  sta{day} = cell(1,num_units);
  for unit=1:num_units
    cfg.spikechannel = spike_labs{unit};
    sta_out = ft_spiketriggeredaverage(cfg, ft_data);
    time = sta_out.time;
    sta{day}{unit} = squeeze(mean(sta_out.trial,2))'; % average over electrode
    if find(size(sta{day}{unit}) == 1)
      sta{day}{unit} = sta{day}{unit}';
    end
  end
  sta{day} = [sta{day}{:}];
end
sta = [sta{:}];
sta = sta(:,~any(isnan(sta),1));
std_err = 1.96 * std(sta,0,2) ./ sqrt(size(sta,2));
sta = mean(sta, 2);
end
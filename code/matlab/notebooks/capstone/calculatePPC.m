function [ppcs, freqs] = calculatePPC(lfps, spikes, electrode_mappings, trial_beg, Fs)
num_days = numel(lfps);
if ~iscell(electrode_mappings)
  electrode_mappings = {electrode_mappings};
end

spec_cfg = [];
spec_cfg.method = 'mtmfft';
spec_cfg.foilim = [5, 100];
spec_cfg.timwin = [-0.05, 0.05]; % around each spike
spec_cfg.taper = 'hanning'; % or dpss for multitaper
spec_cfg.feedback = 'no';

stat_cfg = [];
method = 'ppc0'; % or ppc1 or ppc2
stat_cfg.method = method; 
stat_cfg.timwin = 'all'; % don't make a "spectrogram" of ppc values
stat_cfg.avgoverchan = 'unweighted';
stat_cfg.latency = 'prestim'; % only calculate PPC before shape onset

ppcs = cell(1,num_days);
freqs = cell(1,num_days);
for day=1:num_days
  fprintf("Day %d / %d\n", day, num_days);
  [ft_data, lfp_labels, spike_labels] = toFieldtrip(lfps{day}, spikes{day}, trial_beg, Fs);
  ppcs{day} = cell(1,numel(spike_labels));
  freqs{day} = cell(1,numel(spike_labels));
  for unit_num=1:numel(spike_labels)
    spec_cfg.spikechannel = spike_labels{unit_num};
    spec_out = ft_spiketriggeredspectrum(spec_cfg, ft_data);
    
    exclude_idx = getElectrodeIdxFromUnitNum(electrode_mappings{day}, unit_num);
    exclude_chan = electrode_mappings{day}{exclude_idx,1};
    chans = true(1,numel(lfp_labels));
    chans(exclude_chan) = false;
    
    stat_cfg.spikechannel = spike_labels{unit_num};
    stat_cfg.channel = lfp_labels(chans);
    stat_out = ft_spiketriggeredspectrum_stat(stat_cfg, spec_out);
    if isnan(stat_out.(method))
      ppcs{day}{unit_num} = []; % sometimes it just returns a single day, so this prevents it from breaking
      freqs{day}{unit_num} = [];
    else
      ppcs{day}{unit_num} = stat_out.(method)';
      freqs{day}{unit_num} = stat_out.freq';
    end
    freqs{day}{unit_num} = stat_out.freq';
  end
  ppcs{day} = [ppcs{day}{:}];
  freqs{day} = [freqs{day}{:}];
end
end
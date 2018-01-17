function plt = plotSample(dh, start_time, num_seconds, unit_num)
%PLOTSAMPLE For use in data exploration and verification
% Plots the mean of the LFPs for a number of seconds (num_seconds) after
% the start time (start_time). Alongside the mean LFP, it plots the spikes
% from a specified unit and the trial events in different colors.
plt = figure('Visible', 'off');
plot(mean(dh.lfps, 2));
xlim([start_time start_time+num_seconds*dh.lfp_sample_freq]);
ylim([-4000 4000]);
times = string(duration(0,0,0,1:dh.lfp_sample_freq:size(dh.lfps,1), 'Format', 'mm:ss'));
set(gca, 'XTick', 1:dh.lfp_sample_freq:size(dh.lfps,1), 'XTickLabel', times);
hold on

for i=1:numel(dh.noise), plot([dh.noise(i) dh.noise(i)], [-2000 2000], 'Color', 'red', 'LineWidth', 1); end
for i=1:numel(dh.fixate), plot([dh.fixate(i) dh.fixate(i)], [-2000 2000], 'Color', 'green', 'LineWidth', 1); end
for i=1:numel(dh.shape), plot([dh.shape(i) dh.shape(i)], [-2000 2000], 'Color', 'yellow', 'LineWidth', 1); end
for i=1:numel(dh.saccade), plot([dh.saccade(i) dh.saccade(i)], [-2000 2000], 'Color', 'magenta', 'LineWidth', 1); end

if nargin > 3
  spikes = dh.spikes{unit_num};
  spikes = spikes(and(spikes >= start_time, spikes < start_time+num_seconds*dh.lfp_sample_freq));
  for i=1:numel(spikes), plot([spikes(i) spikes(i)], [-2000 2000], 'Color', 'black', 'LineWidth', 0.1); end
end
end
function plts = plotCombined(dh, between, unit_number)
plts = gobjects(1,6);
window_sizes = [50 100 200 501];
ticks = round(linspace(1,diff(between),10));
t = 1:diff(between);

mean_lfp = mean(dh.lfps,2);
mean_lfp = mean_lfp(between(1):between(2)-1);
z = hilbert(mean_lfp);
envelope = abs(z);
phase = unwrap(angle(z));
freq = diff(phase)*dh.lfp_sample_freq/(2*pi);

spikes = dh.spikes{unit_number};
spikes = spikes(and(spikes >= between(1), spikes < between(2)));
spike_signal = zeros(1, numel(mean_lfp));
spike_signal(spikes-between(1)+1) = 1;

plts(1) = figure('Visible', 'off');
plot(t(2:end), freq);
title('LFP Instantaneous Frequency');
setTicks(ticks, between(1));

plts(2) = figure('Visible', 'off');
plot(t, envelope);
title('LFP Envelope');
setTicks(ticks, between(1));

plts(3) = figure('Visible', 'off');
default_color_order = get(gcf, 'DefaultAxesColorOrder');
set(gcf, 'DefaultAxesColorOrder', parula(numel(window_sizes)));
hold on
legend show
for size=window_sizes
  window = gausswin(size);
  smoothed = conv(window, spike_signal);
  plot(t(ceil(size/2):end-floor(size/2)), smoothed(size:end-size+1)/sum(window), 'DisplayName', sprintf('%d ms', size));
end
title(sprintf('Smoothed Spike Rates by Binwidth for Unit %d', unit_number));
setTicks(ticks, between(1));
hold off
set(gcf, 'DefaultAxesColorOrder', default_color_order);

plts(4) = figure('Visible', 'off');
hold on
for i=1:numel(dh.spikes)
  plot_spikes = dh.spikes{i}(and(dh.spikes{i} >= between(1), dh.spikes{i} < between(2)));
  plotLinesAt(plot_spikes, between(1), i-1);
end
title('Spike Trains');
ylim([0 numel(dh.spikes)]);
setTicks(ticks, between(1));
hold off

plts(5) = figure('Visible', 'off');
spectrogram(mean_lfp, hann(128), 64, 512, dh.lfp_sample_freq, 'yaxis');
ylim([0 100]);

plts(6) = figure('Visible', 'off');
plot(t, mean_lfp, 'Color', 'black');
hold on
fixate = dh.fixate(and(dh.fixate >= between(1), dh.fixate < between(2)));
noise = dh.noise(and(dh.noise >= between(1), dh.noise < between(2)));
shape = dh.shape(and(dh.shape >= between(1), dh.shape < between(2)));
saccade = dh.saccade(and(dh.saccade >= between(1), dh.saccade < between(2)));
plot([fixate; fixate]-between(1)+1, [-1000;1000], 'Color', 'red');
plot([noise; noise]-between(1)+1, [-1000;1000], 'Color', 'green');
plot([shape; shape]-between(1)+1, [-1000;1000], 'Color', 'blue');
plot([saccade; saccade]-between(1)+1, [-1000;1000], 'Color', 'magenta');
setTicks(ticks, between(1));
hold off

  function setTicks(ticks, start)
    times = string(duration(0, 0, 0, ticks+start-1, 'Format', 'mm:ss'));
    set(gca, 'XTick', ticks, 'XTickLabel', times);
  end
  function plotLinesAt(times, start, bottom)
    plot([times; times]-start+1, [0; 1]+bottom, 'Color', 'black', 'LineWidth', 0.1)
  end
end
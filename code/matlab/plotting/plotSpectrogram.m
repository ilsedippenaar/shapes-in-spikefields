function plt = plotSpectrogram(S,t,f,take_log)
if nargin <= 3 || take_log
  take_log = true;
  S = 10*log10(S);
end
plt = figure('Visible', 'off');
imagesc(flipud(S'));
bar = colorbar;
colormap hot;
if take_log
  ylabel(bar, 'Power (dB)');
end

times = linspace(t(1),t(end),10);
time_ticks = interp1(times, linspace(1,numel(t),10), times);
time_labels = arrayfun(@(x) sprintf('%.0f',x), times, 'UniformOutput', false);
set(gca, 'XTick', time_ticks, 'XTickLabel', time_labels);
xlabel('Time (s)');

freqs = linspace(f(end),f(1),10);
freq_ticks = interp1(freqs, linspace(1,numel(f),10), freqs);
freq_labels = arrayfun(@(x) sprintf('%.1f',x), freqs, 'UniformOutput', false);
set(gca, 'YTick', freq_ticks, 'YTickLabel', freq_labels);
ylabel('Frequency (Hz)');
end
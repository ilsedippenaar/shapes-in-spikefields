function plt = plotSpectrogram(S,t,f)
plt = figure('Visible', 'off');
imagesc(flipud(10*log10(S)'));
bar = colorbar;
colormap hot;
ylabel(bar, 'Power (dB)');

times = linspace(t(1),t(end),10);
time_ticks = interp1(times, linspace(1,numel(t),10), times);
time_labels = arrayfun(@(x) sprintf('%.2f',x), times, 'UniformOutput', false);
set(gca, 'XTick', time_ticks, 'XTickLabel', time_labels);
xlabel('Time (s)');

freqs = linspace(f(end),f(1),10);
freq_ticks = interp1(freqs, linspace(1,numel(f),10), freqs);
freq_labels = arrayfun(@(x) sprintf('%.1f',x), freqs, 'UniformOutput', false);
set(gca, 'YTick', freq_ticks, 'YTickLabel', freq_labels);
ylabel('Frequency (Hz)');
end
function plt = electrodeHeatmap(mat, plot_ax)
assert(all(size(mat) == [10 10]))
range = [min(min(mat)) max(max(mat))];
color_vals = linspace(range(1), range(2), 4);
hot_colors = [0 0 0; 1 0 0; 1 1 0; 1 1 1];

plt = figure('Visible', 'off');
if nargin == 1
  plot_ax = axes('Position', [0.05 0.05 0.75 0.8]);
end
set(plot_ax, 'XTick', [], 'YTick' ,[], 'Visible', 'off');
for i=1:10
  for j=1:10
    if ~isnan(mat(i,j))
      color = interp1(color_vals, hot_colors, mat(i,j));
    else
      color = [0 1 1];
    end
    rectangle(plot_ax, 'Position', [j-1, 10-i, 1, 1], ...
      'FaceColor', color, 'EdgeColor', 'none');
  end
end
ticks = interp1([color_vals(1) color_vals(end)], [0 1], color_vals);
tick_labels = arrayfun(@(x) sprintf('%.1f',x), color_vals, 'UniformOutput', false);
if nargin == 1 % colorbar position is annoying with subplots
  colorbar(plot_ax, 'Ticks', ticks, 'TickLabels', tick_labels);
  colormap(plot_ax, interp1(color_vals, hot_colors, linspace(color_vals(1), color_vals(end), 1000)));
end
end
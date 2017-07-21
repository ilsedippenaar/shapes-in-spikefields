function plt = electrodeHeatmap(mat)
assert(all(size(mat) == [10 10]))
range = [min(min(mat)) max(max(mat))];
color_vals = linspace(range(1), range(2), 4);
hot_colors = [0 0 0; 1 0 0; 1 1 0; 1 1 1];

square_size = 100;
plt = figure('Visible', 'off');
set(gca, 'Position', [0 0 1 1], 'Visible', 'off');
plot_ax = axes('Position', [0.05 0.05 0.75 0.8]);
set(plot_ax, 'XTick', [], 'YTick' ,[]);
for i=1:10
  for j=1:10
    if ~isnan(mat(i,j))
      color = interp1(color_vals, hot_colors, mat(i,j));
    else
      color = [0 1 1];
    end
    rectangle(plot_ax, 'Position', [j-1, 10-i, 1, 1]*square_size, ...
      'FaceColor', color, 'EdgeColor', 'none');
  end
end
ticks = interp1([color_vals(1) color_vals(end)], [0 1], color_vals);
tick_labels = arrayfun(@(x) sprintf('%.1f',x), color_vals, 'UniformOutput', false);
colorbar('Position', [0.825, 0.1, 0.05, 0.7], ...
  'Ticks', ticks, 'TickLabels', tick_labels);
colormap(interp1(color_vals, hot_colors, linspace(color_vals(1), color_vals(end), 1000)));
end
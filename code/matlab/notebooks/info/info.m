%% Calculate information surfaces
bws = 25:25:250;
delays = 0:25:500;
data_types = {'lfp', 'spike'};
info_types = {'sensory', 'choice'};
combs = cellfun(@(c) [c{:}], cartesianProd(data_types, info_types), 'UniformOutput', false);
combs{end+1} = {[], 'behavioral'};
info_mats = cell(1,numel(combs));
for i=1:numel(combs)
  fprintf('%s %s\n', combs{i}{1}, combs{i}{2});
  info_mats{i} = funcOverVecs(@(bw,d) ...
    calculateInfo(dhs, combs{i}{1}, combs{i}{2}, bw, d, 'average', true), bws, delays);
end
%% Plotting
for i=1:numel(combs)
  plt = plotInfoSurface(info_mats{i}, combs{i}{1}, combs{i}{2}, bws, delays);
  if isempty(combs{i}{1})
    saveFigures(plt, fullfile(plot_save_dir, 'info', sprintf('%s.png', combs{i}{2})));
  else
    saveFigures(plt, fullfile(plot_save_dir, 'info', sprintf('%s_%s.png', combs{i}{:})));
  end
end
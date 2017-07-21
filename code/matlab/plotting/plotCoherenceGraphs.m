function plts = plotCoherenceGraphs(dh, spec, start_time, end_time, varargin)
p = inputParser;
p.addParameter('freq_bin', [0 100], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('window_width', 2^(nextpow2(dh.lfp_sample_freq)));
p.addParameter('num_fft', 2^(nextpow2(dh.lfp_sample_freq)));
p.parse(varargin{:});
args = p.Results;

lfps = dh.getLfpSlices(spec, start_time, end_time, 'data_type', 'single');
plts = gobjects(1,dh.num_lfp_electrodes);
for i=1:dh.num_lfp_electrodes
  r = int32(rem(numel(lfps{i}), args.window_width)); % this is all to prevent what is presumably a bug in mscohere
  plts(i) = figure('Visible', 'off', 'Units', 'points', 'Position', [0 0 2000 2000], 'Color', 'white');
  set(gca, 'Visible', 'off');
  for j=1:dh.num_lfp_electrodes
    str = sprintf('%d of %d', j, dh.num_lfp_electrodes);
    fprintf('%s',str);
    [cxy,f] = mscohere(lfps{i}(1:end-r,:), lfps{j}(1:end-r, :), hann(args.window_width),...
      args.window_width/2, args.num_fft, dh.lfp_sample_freq);
    electrode_idx = find(cellfun(@(c) ~isempty(c) && c == j, {dh.electrode_mapping{:,1}}));
    if electrode_idx < 9
      plt_idx = electrode_idx + 1;
    elseif electrode_idx < 89
      plt_idx = electrode_idx + 2;
    else
      plt_idx = electrode_idx + 3;
    end
    ax = subplot(10,10,plt_idx);
    freq_idxs = and(f >= args.freq_bin(1), f < args.freq_bin(2));
    plotMeanAndStds(cxy(freq_idxs,:), 'x', f(freq_idxs));
    ax.Color = getColorFromCoherence(mean(mean(cxy(freq_idxs,:),2)));
    fprintf(sprintf('%s', strjoin(repmat({'\b'},1,numel(str))),'')); % why matlab whyyyyy
  end
  pos = get(subplot(10,10,100, 'Visible','off'), 'Position');
  ticks = [0 0.5 0.8 0.9 0.95 1];
  colorbar(gca, 'Position', [pos(1)+pos(3)*1.1, 0.2, pos(3)/3, 0.6], ...
    'Ticks', ticks, 'TickLabels', cellstr(string(ticks)));
  colormap(gca, getColorFromCoherence(0:0.001:1));
  fprintf('Finished electrode %d, elapsed time = %.2f\n', i, toc);
end

  function col = getColorFromCoherence(coh)
    colors = [0 0 0; ...
              1 0 0; ...
              1 1 0; ...
              1 1 1];
    col = interp1([0 0.9 0.95 1], colors, coh);
  end
end
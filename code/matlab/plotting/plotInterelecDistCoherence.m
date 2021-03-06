function plt = plotInterelecDistCoherence(cohs, freqs, num_in_bins, dist_bins, varargin)
p = inputParser;
p.addParameter('freq_bin', [min(freqs) max(freqs)+1]);
p.addParameter('std_err', []);
p.parse(varargin{:});
args = p.Results;

plt = figure('Visible', 'off');
n = size(cohs, 2);
hold on
set(gca, 'ColorOrder', flipud(jet(n)));
idxs = and(freqs >= args.freq_bin(1), freqs < args.freq_bin(2));
for i=1:n
  name = sprintf('%.0f - %.0f $\\mu m$, n=%d', dist_bins(1,i), dist_bins(2,i), num_in_bins(i));
  plot(freqs(idxs), cohs(idxs,i), 'LineWidth', 2, 'DisplayName', name);
  if args.std_err
    errorbar(freqs(idxs), cohs(idxs,i), args.std_err(idxs,i));
  end
end
l = legend('show');
set(l, 'Interpreter', 'latex');
title('Coherences by Interelectrode Distance');
end
function plt = plotInterelecDistCoherence(cohs, freqs, num_in_bins, dist_bins, varargin)
p = inputParser;
p.addParameter('freq_bin', [min(freqs) max(freqs)+1]);
p.parse(varargin{:});
args = p.Results;

plt = figure('Visible', 'off');
n = size(cohs, 2);
hold on
set(gca, 'ColorOrder', flipud(jet(n)));
l = legend('show');
set(l, 'Interpreter', 'latex');
idxs = and(freqs >= args.freq_bin(1), freqs < args.freq_bin(2));
for i=1:n
  name = sprintf('%.0f - %.0f $\\mu m$, n=%d', dist_bins(1,i), dist_bins(2,i), num_in_bins(i));
  plot(freqs(idxs), cohs(idxs,i), 'LineWidth', 2, 'DisplayName', name);
end
title('Coherences by Interelectrode Distance');
end
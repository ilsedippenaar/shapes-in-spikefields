function [plts, all_cohs, freqs]=plotCoherenceVariance(cohs, freqs)
all_cohs = zeros(size(cohs, 1), nchoosek(size(cohs,2), 2));
idx = 1;
for i=1:size(cohs,2)
  for j=i+1:size(cohs,3)
    all_cohs(:,idx) = cohs(:,i,j);
    idx = idx + 1;
  end
end

plts = gobjects(1,2);
plts(1) = plotMeanAndStds(all_cohs, 'x', freqs);
title('Mean of all coherences');
xlabel('Frequency (Hz)');
ylabel('Coherence');

plts(2) = figure('Visible', 'off');
plot(freqs, std(all_cohs,0,2) / sqrt(size(all_cohs,2)-1));
title(sprintf('SEM of coherence, N=%d', size(all_cohs,2)));
xlabel('Frequency (Hz)');
ylabel('SEM');
end
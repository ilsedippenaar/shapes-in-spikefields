function varargout=plotCoherenceVariance(cohs, freqs)
all_cohs = zeros(size(cohs, 1), nchoosek(size(cohs,2), 2));
idx = 1;
for i=1:size(cohs,2)
  for j=i+1:size(cohs,3)
    all_cohs(:,idx) = cohs(:,i,j);
    idx = idx + 1;
  end
end

if nargout ~= 0
  fig = figure('Visible', 'off');
else
  figure;
end
plotMeanAndStds(all_cohs, 'x', freqs);
title('Mean of all coherences');
xlabel('Frequency (Hz)');
ylabel('Coherence');

if nargout ~= 0
  fig_stds = figure('Visible', 'off');
else
  figure;
end
plot(freqs, std(all_cohs,0,2) / sqrt(size(all_cohs,2)-1));
title('SEM of coherence');
xlabel('Frequency (Hz)');
ylabel(sprintf('SEM, N=%d', size(all_cohs,2)));

if nargout >= 1
  varargout{1} = [fig fig_stds];
  if nargout >= 2
    varargout{2} = all_cohs;
    if nargout == 3
      varargout{3} = freqs;
    end
  end
end
end
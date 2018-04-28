function phases = getHHTPhases(x, freq_bins, phase_at, freq_between, Fs)
if nargin < 5
  Fs = 1000;
end
window_size = freq_between(2) - freq_between(1) + 1;
idx_to_freq = (Fs/2) / numel(pwelch(x(freq_between(1):freq_between(2))), hann(window_size));
imfs = emd(x);
imf_freqs = cellfun(@(c) argmax(pwelch(c(freq_between(1):freq_between(2)), hann(window_size))) * idx_to_freq, imfs); % get frequency with greatest power in each IMF
phases = cell(1,numel(freq_bins)-1);
for i=1:numel(freq_bins)-1
  freq_imfs = imfs(and(imf_freqs >= freq_bins(i), imf_freqs < freq_bins(i+1)));
  phases{i} = zeros(1,numel(freq_imfs));
  for j=1:numel(freq_imfs)
    h = hilbert(freq_imfs{j});
    if phase_at == 'end'
      phases{i}(j) = h(end);
    else
      phases{i}(j) = h(phase_at);
    end
  end
end
end
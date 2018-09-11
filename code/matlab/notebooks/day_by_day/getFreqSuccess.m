function [x, y, std_err] = getFreqSuccess(tbl, label_name, analysis_idxs, freq_cutoff, Fs)
freqs = getFourierFreqs(analysis_idxs, Fs);
success_counts = countcats(categorical(tbl(tbl.success,:).(label_name)));
failure_counts = countcats(categorical(tbl(~tbl.success,:).(label_name)));

success_rate = success_counts ./ (success_counts+failure_counts);
x = freqs(freqs < freq_cutoff);
y = success_rate(:, freqs < freq_cutoff)';
% take binomial (bernoulli) error for success rates
n = success_counts + failure_counts;
std_err = sqrt(success_rate .* (1-success_rate) ./ n) * 1.96;
std_err = std_err(:, freqs < freq_cutoff)';
end
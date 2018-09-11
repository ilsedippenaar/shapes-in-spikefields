function [S,t,f] = calculateHilbertSpectrum(X, Fs)
% X is a vector
% this function accounts for possible negative frequencies
X = (X-mean(X))/var(X);
[Q, D] = getSSAComponents(X, Fs);
D = diag(D)/sum(D);

z = hilbert(Q*D);

freq = diff(unwrap(angle(z))) * Fs / (2*pi); % instantaneous frequency
neg_freq_bins = getFourierFreqs(1:numel(X),-Fs);
neg_freq_bins = neg_freq_bins(end:-1:2); % exclude 0
freq_bins = [neg_freq_bins, getFourierFreqs(1:numel(X),Fs)];
freq_labs = discretize(freq, freq_bins);

mags = movmean(abs(z), 2, 'Endpoints', 'discard'); % this is necessary to line magnitudes up with instaneous frequency

subs = [repmat((1:size(z,1)-1)', size(z,2), 1), ...
        freq_labs(:)];
S = accumarray(subs, mags(:), [size(mags,1),size(freq_bins,2)-1]); % time x [-freq,freq] (should be almost square)
t = (1:numel(X)-1)'/Fs;
f = movmean(freq_bins,2,'EndPoints', 'discard')';
end
function [Q, D, f] = getSSAComponents(v, Fs)
v = (v - mean(v)) / std(v);
[Q,D] = ssacom(v, round(Fs/10)); % HACK: can detect a 10 Hz signal
S = abs(fft(Q, [], 1));
[~, idxs] = max(S);
f = (idxs-1) * Fs / numel(v);
%c = Q(:,f < 100)*diag(D(f<100))/sum(D); % less than 100 Hz
%f = f(f<100);
%D = D(f<100) / sum(D);
end
Fs = 1000;
T = 3;
N = T*Fs;
M = 3;

t = linspace(0,T,N);
f = [3; 12; 25];
ft = repelem(f, 1, N) + diag(f/100) * randn(M,N) ./ repelem(t+1, M, 1);
a = [2; 3; 5];
at = repelem(a, 1, N) + diag(a/100) * randn(M,N);
p = [0; pi/2; 5*pi/6];
X = sum(at .* cos(2*pi*ft.*repelem(t,M,1) + repelem(p, 1, N)))' + randn(N,1);
X = (X-mean(X))/var(X);
[Q, D, freqs] = getSSAComponents(X, Fs);
disp('Freqs and percent contribution');
disp([freqs(1:10);(D(1:10)'/sum(D))])
D = diag(D)/sum(D);
% idxs = 1:M*2;
idxs = 1:size(D,2);

figure
recon = sum(Q(:,idxs)*D(idxs,idxs),2);
plot(t, [X,recon])
xlim([0,0.5])
title(sprintf('SSA Reconstuction with %d components', M));
legend('Original signal', 'Reconstructed');

figure
pwelch(recon,[],[],[],Fs)
title('PSD of reconstructed signal');
[pks,locs] = baseFindpeaks(log(pwelch(recon, hann(2^nextpow2(N/4)))), ...
  'NPeaks', M, 'SortStr', 'descend', 'MinPeakProminence', 1);
disp('Ground truth and estimated frequencies in reconstruction');
disp([f'; interp1([0,256], [0,Fs/2], locs')])
% frequency content as expected

bpows_X = arrayfun(@(x) bandpower(X,Fs,[x-2,x+2]),f);
bpows_recon = arrayfun(@(x) bandpower(recon,Fs,[x-2,x+2]),f);
disp('Ground truth magnitudes, bandpower of orignal signal in frequency bins, bandpower of reconstructed');
disp([a';bpows_X'/max(bpows_X)*max(a);bpows_recon'/max(bpows_recon)*max(a)])
% magnitude of pwelch doesn't match original contributions, even for the
% original signal
% note: plain old fft works here, but more difficult to find peaks

i = find(t>=0.25,1);
z = hilbert(Q(:,idxs)*D(idxs,idxs));
disp('Ground truth phase angle, phase angle of components');
disp(mod(0.25*f'*2*pi+p',2*pi)/(2*pi))
disp(mod(angle(z(i,:)),2*pi)/(2*pi))
% close enough phase angles to expected

figure
plot(abs(z(:,idxs)))
title('Envelope of reconstructed components');
mags = mean(abs(z(:,idxs)));
disp('Ground truth magnitudes, average envelope magnitude of reconstructed components');
disp(a')
disp(mags/max(mags)*max(a))
% not even close on average envelope - this is probably an entirely wrong
% technique to approximate frequency content

% Hilbert Spectrum
[S,t_spec,f_spec] = calculateHilbertSpectrum(X,Fs); % t x f

figure
imagesc(10*log10(flipud(S')))
title('Hilbert Spectrum');
colorbar
colormap(flipud(hot(50)));

set(gca,'XTick', 1:300:numel(t_spec), 'XTickLabel', num2str(t_spec(1:300:end),2));

f_spec = flipud(f_spec);
set(gca, 'YTick', 1:50:numel(f_spec), 'YTickLabel', f_spec(1:50:end))
ylim([find(f_spec<40,1), find(f_spec<0,1)]) % imagesc counts down, so y=0 is at the top

figure
plot(flipud(f_spec), 10*log10(mean(S)'))
title('Marginal Hilbert Spectrum');
xlim([0,Fs/2])
xlabel('Time (s)');
ylabel('Power (AU)');
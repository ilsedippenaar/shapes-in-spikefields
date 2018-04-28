function phases = getLfpPhases(lfps, freq_bins, phase_at, freq_between, Fs, method)
% lfps is cell array with days -> electrodes -> time x trials int16
if nargin < 6
  method = 'hht';
end
phases = cell(1,numel(lfps));
for i=1:numel(lfps) % day
  day_lfps = lfps{i}(~cellfun(@isempty, lfps{i})); % eliminate empty electrodes
  if isempty(day_lfps)
    phases{i} = {};
    continue
  end
  day_lfps = reshape([lfps{i}{:}], [size(day_lfps{1}), numel(day_lfps)]); % time x trials x electrode
  day_lfps = mean(day_lfps,3); % time x trials
  phases{i} = cell(1,size(day_lfps,2));
  for j=1:size(day_lfps,2) % for each trial
    if strcmp(method, 'hht')
      phases{i}{j} = getHHTPhases(day_lfps(:,j), freq_bins, phase_at, freq_between, Fs);
    elseif strcmp(method, 'fft')
      N = diff(freq_between) + 1;
      f = Fs / N  * (-N/2:(N/2)-1);
      out = fftshift(fft(day_lfps(freq_between(1):freq_between(2),j)));
      out = out(f>=0);
      f = f(f>=0);
      phases{i}{j} = cell(1,numel(freq_bins)-1);
      for k=1:numel(freq_bins)-1
        phases{i}{j}{k} = out(and(f >= freq_bins(k), f < freq_bins(k+1)));
      end
    end
  end
  phases{i} = combineCellArrays('complex', phases{i}{:});
end
phases = combineCellArrays('complex', phases{:});
end
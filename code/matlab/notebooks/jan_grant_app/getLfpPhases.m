function phases = getLfpPhases(lfps, freq_bins, phase_at, freq_between, Fs)
% lfps is cell array with days -> electrodes -> time x trials int16
phases = cell(1,numel(lfps));
for i=1:numel(lfps) % day
  day_lfps = lfps{i}(~cellfun(@isempty, lfps{i})); % eliminate empty electrodes
  day_lfps = reshape([lfps{i}{:}], [size(day_lfps{1}), numel(day_lfps)]); % time x trials x electrode
  day_lfps = mean(day_lfps,3); % time x trials
  phases{i} = cell(1,size(day_lfps,2));
  for j=1:size(day_lfps,2) % for each trial
    phases{i}{j} = getHHTPhases(day_lfps(:,j), freq_bins, phase_at, freq_between, Fs);
  end
  phases{i} = combineCellArrays('complex', phases{i}{:});
end
phases = combineCellArrays('complex', phases{:});
end
function plt = plotPsdHeatmap(dh, psds, freqs, freq_bin)
psd_map = nan(10);
idxs = and(freqs >= freq_bin(1), freqs < freq_bin(2));
for i=1:size(psds,2)
  if i == 55 % this is a bad electrode
    continue
  end
  pos = getLocFromElectrodeIdx(getElectrodeIdxFromLfpIdx(dh, i));
  psd_map(pos(1), pos(2)) = mean(psds(idxs,i));
end
plt = electrodeHeatmap(psd_map);
title(sprintf('%d to %d Hz', freq_bin(1), freq_bin(2)));
end
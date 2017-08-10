function plt = plotPsdHeatmap(psds, freqs, freq_bin, electrode_mapping)
psd_map = nan(10);
idxs = and(freqs >= freq_bin(1), freqs < freq_bin(2));
for i=1:size(psds,2)
  pos = getLocFromElectrodeIdx(getElectrodeIdxFromLfpIdx(electrode_mapping, i));
  psd_map(pos(1), pos(2)) = mean(psds(idxs,i));
end
plt = electrodeHeatmap(psd_map);
title(sprintf('%d to %d Hz', freq_bin(1), freq_bin(2)));
end
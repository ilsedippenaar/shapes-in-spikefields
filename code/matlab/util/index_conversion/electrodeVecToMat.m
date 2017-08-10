function mat = electrodeVecToMat(electrode_mapping, vec)
mat = nan(10);
for i=1:numel(vec)
  pos = getLocFromElectrodeIdx(getElectrodeIdxFromLfpIdx(electrode_mapping, i));
  mat(pos(1), pos(2)) = vec(i);
end
end
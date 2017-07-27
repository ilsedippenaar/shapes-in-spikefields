function mat = electrodeVecToMat(dh, vec)
mat = nan(10);
for i=1:numel(vec)
  pos = getLocFromElectrodeIdx(getElectrodeIdxFromLfpIdx(dh, i));
  mat(pos(1), pos(2)) = vec(i);
end
end
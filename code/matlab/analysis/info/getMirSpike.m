function mir = getMirSpike(d, bw)
if isempty(d)
  mir = [];
  return
end
mir = mutInfo(~cellfun(@isempty, d.data), grp2idx(d.label)) * 1000 / bw;
end
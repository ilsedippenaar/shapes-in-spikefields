function out = expandCellArray(electrode_mapping, c)
% c is a cell array of electrode data
out = cell(1,96);
for i=1:numel(c)
  elec_idx = getElectrodeIdxFromLfpIdx(electrode_mapping, i);
  out{elec_idx} = c{i};
end
end
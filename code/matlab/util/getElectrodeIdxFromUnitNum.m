function elec_idx = getElectrodeIdxFromUnitNum(dh, unit_num)
elec_idx = find(cellfun(@(c) ~isempty(c) && any(c==unit_num), dh.electrode_mapping(:,2)));
end
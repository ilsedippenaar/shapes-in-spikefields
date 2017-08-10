function elec_idx = getElectrodeIdxFromUnitNum(electrode_mapping, unit_num)
elec_idx = find(cellfun(@(c) ~isempty(c) && any(c==unit_num), electrode_mapping(:,2)));
end
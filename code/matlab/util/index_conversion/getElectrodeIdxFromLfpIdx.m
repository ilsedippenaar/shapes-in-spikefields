function elec_idx = getElectrodeIdxFromLfpIdx(electrode_mapping, lfp_idx)
elec_idx = find(cellfun(@(c) ~isempty(c) && c==lfp_idx, electrode_mapping(:,1)));
end
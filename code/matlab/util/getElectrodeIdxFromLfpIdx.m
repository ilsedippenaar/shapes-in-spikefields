function elec_idx = getElectrodeIdxFromLfpIdx(dh, lfp_idx)
elec_idx = find(cellfun(@(c) ~isempty(c) && c==lfp_idx, dh.electrode_mapping(:,1)));
end
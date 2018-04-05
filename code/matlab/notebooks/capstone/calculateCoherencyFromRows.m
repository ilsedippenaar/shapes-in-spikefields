function cohs = calculateCoherencyFromRows(lfps, electrode_mapping, mt_params, analysis_window_idxs)
% horrid series of manipulations to get lfps from trial -> time x
% electrodes to day -> electrode -> time x trial
lfps = cellfun(@(x) mat2cell(single(x), size(x,1), ones(1,size(x,2))), lfps, 'UniformOutput', false);
lfps = {cellfun(@(c) [c{:}], zipCell(lfps{:}), 'UniformOutput', false)};

cohs = calculateCoherency(lfps, electrode_mapping(1), mt_params, analysis_window_idxs);
end
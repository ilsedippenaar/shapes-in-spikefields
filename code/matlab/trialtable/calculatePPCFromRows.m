function [ppcs, freqs] = calculatePPCFromRows(lfps, spikes, electrode_mapping, trial_beg, Fs)
% horrid series of manipulations to get lfps from trial -> time x
% electrodes to day -> electrode -> time x trial
lfps = cellfun(@(x) mat2cell(single(x), size(x,1), ones(1,size(x,2))), lfps, 'UniformOutput', false);
lfps = {cellfun(@(c) [c{:}], zipCell(lfps{:}), 'UniformOutput', false)};

spikes = {zipCell(spikes{:})};
[ppcs, freqs] = calculatePPC(lfps, spikes, electrode_mapping(1), trial_beg, Fs);
end
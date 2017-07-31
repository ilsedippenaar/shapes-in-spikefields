function obj = clean(obj, dist_mat, varargin)
p = inputParser;
p.addParameter('cosine_dist_cutoff', 0.5);
p.parse(varargin{:});
args = p.Results;

valid_idxs = find(mean(dist_mat,1) >= args.cosine_dist_cutoff);
valid_elec_idxs = arrayfun(@(x) getElectrodeIdxFromLfpIdx(obj, x), valid_idxs);
valid_spike_idxs = cell2mat(arrayfun(@(x) obj.electrode_mapping{x,2}, valid_elec_idxs, 'UniformOutput', false));
valid_spike_elec_idxs = arrayfun(@(x) getElectrodeIdxFromUnitNum(obj, x), valid_spike_idxs);

fprintf('Eliminating %d electrodes\n', obj.num_lfp_electrodes-numel(valid_idxs));

obj.lfps = obj.lfps(:,valid_idxs);
obj.spikes = obj.spikes(valid_spike_idxs);
obj.unit_snr = obj.unit_snr(valid_spike_idxs);

obj.num_lfp_electrodes = size(obj.lfps,2);
obj.num_units = numel(obj.spikes);

electrode_mapping = cell(size(obj.electrode_mapping));
for i=1:numel(valid_elec_idxs)
  electrode_mapping{valid_elec_idxs(i),1} = i;
end
for i=1:numel(valid_spike_elec_idxs)
  electrode_mapping{valid_spike_elec_idxs(i),2} = [electrode_mapping{valid_spike_elec_idxs(i),2}, i];
end
obj.electrode_mapping = electrode_mapping;
end
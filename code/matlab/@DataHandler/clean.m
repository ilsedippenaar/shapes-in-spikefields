function obj = clean(obj, varargin)
p = inputParser;
p.addParameter('cosine_dist_cutoff', 0.5);
p.addParameter('remove_60', true);
p.addParameter('lowpass', 100);
p.parse(varargin{:});
args = p.Results;

dist_mat = calculateAdjMatrix(obj.lfps, 'cosine');
valid_idxs = find(mean(dist_mat,1) >= args.cosine_dist_cutoff);
valid_elec_idxs = arrayfun(@(x) getElectrodeIdxFromLfpIdx(obj.electrode_mapping, x), valid_idxs);
valid_spike_idxs = cell2mat(arrayfun(@(x) obj.electrode_mapping{x,2}, valid_elec_idxs, 'UniformOutput', false));
valid_spike_elec_idxs = arrayfun(@(x) getElectrodeIdxFromUnitNum(obj.electrode_mapping, x), valid_spike_idxs);

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

% if args.remove_60
%   w0 = 60/(obj.lfp_sample_freq/2);
%   [b,a] = iirnotch(w0, w0/35);
%   obj.lfps = filter(b, a, single(obj.lfps));
% end
% if args.lowpass
%   [z,p,k] = butter(8, args.lowpass / (obj.lfp_sample_freq/2), 'low');
%   obj.lfps = sosfilt(zp2sos(z,p,k), single(obj.lfps));
% end
end
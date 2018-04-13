function cohs = cohByDay(lfps, electrode_mapping, analysis_idxs, mt_params, freq_idxs)
% lfps = time x electrode x trial or 1 x trial -> time x electrode
if iscell(lfps)
  lfps = cat(3,lfps{:});
end
if nargin < 4
  T = 0.5;
  W = 6;
  mt_params = [];
  mt_params.tapers = [T*W, 2*T*W-1];
  mt_params.Fs = 1000;
end
if ~isfield(mt_params, 'trialave')
  mt_params.trialave = false;
end
num_trials = size(lfps,3);
lfps = single(lfps);
elec_idxs = arrayfun(@(x) getElectrodeIdxFromLfpIdx(electrode_mapping, x), 1:size(lfps,2));

C = coherencyc(lfps(analysis_idxs,1,1), lfps(analysis_idxs,1,1), mt_params);
if nargin < 5
  freq_idxs = 1:numel(C);
end
if mt_params.trialave
  cohs = nan(numel(freq_idxs), 96, 96);
  for i=1:size(lfps,2)
    elec_i = elec_idxs(i);
    cohs(:,elec_i,elec_i)=1;
    for j=i+1:size(lfps,2)
      elec_j = elec_idxs(j);
      coh = coherencyc(lfps(analysis_idxs,i,:), lfps(analysis_idxs,j,:), mt_params);
      cohs(:,elec_i,elec_j) = coh(freq_idxs);
      cohs(:,elec_j,elec_i) = cohs(:,elec_i,elec_j);
    end
  end
else
  cohs = nan(numel(freq_idxs), num_trials, 96, 96); 
  for i=1:size(lfps,2)
    elec_i = elec_idxs(i);
    cohs(:,:,elec_i,elec_i)=1;
    for j=i+1:size(lfps,2)
      elec_j = elec_idxs(j);
      coh = coherencyc(lfps(analysis_idxs,i,:), lfps(analysis_idxs,j,:), mt_params);
      cohs(:,:,elec_i,elec_j) = coh(freq_idxs,:);
      cohs(:,:,elec_j,elec_i) = cohs(:,:,elec_i,elec_j);
    end
  end
  cohs = permute(cohs, [1,3,4,2]);
end
end
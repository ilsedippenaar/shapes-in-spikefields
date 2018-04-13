function cohs = cohByRow(lfps, electrode_mapping, mt_params)
if nargin < 3
  T = 0.5;
  W = 6;
  mt_params = [];
  mt_params.tapers = [T*W, 2*T*W-1];
  mt_params.Fs = 1000;
end
lfps = single(lfps);
elec_idxs = arrayfun(@(x) getElectrodeIdxFromLfpIdx(electrode_mapping, x), 1:size(lfps,2));
C = coherencyc(lfps(:,1), lfps(:,1), mt_params);
cohs = nan(numel(C), 96, 96);
for i=1:size(lfps,2)
  elec_i = elec_idxs(i);
  for j=i+1:size(lfps,2)
    elec_j = elec_idxs(j);
    cohs(:,elec_i,elec_j) = coherencyc(lfps(:,i), lfps(:,j), mt_params);
    cohs(:,elec_j,elec_i) = cohs(:,elec_i,elec_j);
  end
end
end
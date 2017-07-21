function [inter_cohs, freqs] = calculateInterelecDistCoherence(dh, lfp_slices, dist_bins, cache_dir, conditions, varargin)
[cohs,freqs] = calculateCoherence(lfp_slices, cache_dir, conditions, varargin{:}); % hopefully cached
n = numel(lfp_slices);
inter_cohs = zeros(size(cohs,1), size(dist_bins,2));
for i=1:size(dist_bins,2)
  num = 0;
  for j=1:n
    pos1 = getLocFromElectrodeIdx(getElectrodeIdxFromLfpIdx(dh, j));
    for k=j+1:n
      pos2 = getLocFromElectrodeIdx(getElectrodeIdxFromLfpIdx(dh, k));
      delta = pos1-pos2;
      dist = hypot(delta(1), delta(2)) * 400; % 400 micron between each electrode
      if dist >= dist_bins(1,i) && dist < dist_bins(2,i)
        inter_cohs(:,i) = inter_cohs(:,i) + cohs(:,j,k);
        num = num + 1;
      end
    end
  end
  inter_cohs(:,i) = inter_cohs(:,i) / num;
end
end
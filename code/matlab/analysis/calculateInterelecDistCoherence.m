function [inter_cohs, num_in_bins] = calculateInterelecDistCoherence(cohs, dist_bins, electrode_mapping)
n = size(cohs,2);
inter_cohs = zeros(size(cohs,1), size(dist_bins,2));
num_in_bins = zeros(1,size(dist_bins,2));
for i=1:size(dist_bins,2)
  num = 0;
  for j=1:n
    pos1 = getLocFromElectrodeIdx(getElectrodeIdxFromLfpIdx(electrode_mapping, j));
    for k=j+1:n
      pos2 = getLocFromElectrodeIdx(getElectrodeIdxFromLfpIdx(electrode_mapping, k));
      delta = pos1-pos2;
      dist = hypot(delta(1), delta(2)) * 400; % 400 micron between each electrode
      if dist >= dist_bins(1,i) && dist < dist_bins(2,i)
        inter_cohs(:,i) = inter_cohs(:,i) + cohs(:,j,k);
        num = num + 1;
      end
    end
  end
  inter_cohs(:,i) = inter_cohs(:,i) / num;
  num_in_bins(i) = num;
end
end
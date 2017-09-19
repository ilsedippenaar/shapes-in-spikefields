function [inter_cohs, num_in_bins,binned_cohs] = calculateInterelecDistCoherence(cohs, dist_bins, electrode_mapping)
s = size(cohs);
n = s(end);
inter_cohs = zeros([s(1:end-2),size(dist_bins,2)]);
num_in_bins = zeros(1,size(dist_bins,2));
if nargin > 2
  binned_cohs = cell(size(num_in_bins));
end
for i=1:size(dist_bins,2)
  num = 0;
  for j=1:n
    pos1 = getLocFromElectrodeIdx(getElectrodeIdxFromLfpIdx(electrode_mapping, j));
    for k=j+1:n % exclude identity coherence
      pos2 = getLocFromElectrodeIdx(getElectrodeIdxFromLfpIdx(electrode_mapping, k));
      delta = pos1-pos2;
      dist = hypot(delta(1), delta(2)) * 400; % 400 micron between each electrode
      if dist >= dist_bins(1,i) && dist < dist_bins(2,i)
        if numel(s) == 4
          inter_cohs(:,:,i) = inter_cohs(:,:,i) + cohs(:,:,j,k);
        else
          inter_cohs(:,i) = inter_cohs(:,i) + cohs(:,j,k);
        end
        if nargin > 2
          binned_cohs{i} = [binned_cohs{i}, cohs(:,j,k)];
        end
        num = num + 1;
      end
    end
  end
  if numel(s) == 4
    inter_cohs(:,:,i) = inter_cohs(:,:,i) / num;
  else
    inter_cohs(:,i) = inter_cohs(:,i) / num;
  end
  num_in_bins(i) = num;
end
end
function [sigmas, dists] = getSpatialCohErrors(cohs)
assert(size(cohs,2) == 96 && size(cohs,3) == 96);
distances = cell(96,1);
for i=1:96
  distances{i} = zeros(96-i,1);
  pos1 = getLocFromElectrodeIdx(i);
  for j=i+1:96
    pos2 = getLocFromElectrodeIdx(j);
    distances{i}(j-i) = norm(pos1-pos2) * 0.4; % 0.4 mm between electrodes
  end
end
dists = vertcat(distances{:});
[dists,dist_sort_idxs] = sort(dists);

sigmas = zeros(size(cohs,1),3); % sigma, sigma-c_i, sigma + c_i
for freq_idx=1:size(cohs,1)
  reduced_cohs = cell(96,1);
  for i=1:96
    reduced_cohs{i} = cohs(freq_idx,i+1:end,i).';
  end
  % fit a Gaussian model on even function (i.e. f(-x)=f(x))
  reduced_cohs = vertcat(reduced_cohs{:});
  reduced_cohs = reduced_cohs(dist_sort_idxs); % sort by distance
  d = dists(~isnan(reduced_cohs)); % omitnan in case two electrodes never recorded data together
  reduced_cohs = reduced_cohs(~isnan(reduced_cohs));
  x = [-d; d];
  y = repmat(reduced_cohs, 2, 1);
  % make sure a coefficient is 1 so that model(0)=1 (perfect coherence at distance = 0)
  fitobj = fit(x, y, 'exp(-(x/c)^2)', 'Start', [0, 1]); % mean = 0, std = 1
  conf_ints = confint(fitobj);
  
  sigmas(freq_idx,1) = fitobj.c;
  sigmas(freq_idx,2:3) = conf_ints(:,2);
end
end
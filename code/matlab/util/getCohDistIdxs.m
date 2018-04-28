function idxs = getCohDistIdxs(dist_edges)
if nargin < 1
  dist_edges = linspace(0, hypot(9,9)*400, 6);
end

nbins = numel(dist_edges)-1;
dist_mat = funcOverVecs(@dist, 1:96, 1:96);
dist_labs = discretize(dist_mat, dist_edges);
mask = logical(tril(ones(96)) - eye(96));
idxs = cell(1,nbins);
for i=1:nbins
  idx_mat = and(dist_labs == i, mask);
  [elec_i, elec_j] = ind2sub([96,96], find(idx_mat));
  idxs{i} = [elec_j, elec_i];
end

  function d = dist(x,y)
    d = norm(getLocFromElectrodeIdx(x) - getLocFromElectrodeIdx(y)) * 400;
  end
end
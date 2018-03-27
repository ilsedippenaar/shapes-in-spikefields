function comb = cartesian_prod(varargin)
s = cellfun(@numel, varargin);
% could use ndgrid, but this is easier and more memory efficient
comb = cell(1,prod(s));
idxs = cell(1,numel(varargin));
for i=1:prod(s)
  [idxs{:}] = ind2sub(s, i);
  comb{i} = arrayfun(@(j) varargin{j}(idxs{j}), 1:numel(varargin), 'UniformOutput', false);
end
end
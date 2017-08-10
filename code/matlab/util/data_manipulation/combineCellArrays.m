function out = combineCellArrays(data_type, varargin)
shape = size(varargin{1});
for i=1:numel(varargin)
  assert(iscell(varargin{i}) && all(size(varargin{i}) == shape));
end
out = cell(shape);
for i=1:numel(out)
  c = cell(1,numel(varargin));
  for j=1:numel(varargin)
    c{j} = varargin{j}{i};
  end
  if ~isempty(data_type)
    out{i} = cellArray2mat(c, data_type);
  else
    total = sum(cellfun(@numel, c));
    out{i} = cell(1,total);
    curr_idx = 1;
    for j=1:numel(c)
      n = numel(c{j});
      out{i}(curr_idx:curr_idx+n-1) = c{j};
      curr_idx = curr_idx + n;
    end
  end
end
end
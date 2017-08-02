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
  out{i} = cellArray2mat(c, data_type);
end
end
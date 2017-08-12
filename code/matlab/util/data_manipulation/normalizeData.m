function out = normalizeData(data, varargin)
if ~iscell(data)
  assert(isnumeric(data));
  data = {data};
else
  assert(numel(data) > 0);
end
p = inputParser;
p.addParameter('dim', find(size(data{1}) ~= 1, 1));
p.addParameter('len', []);
p.parse(varargin{:});
args = p.Results;

dim = args.dim;
s = size(data{1});
idxs = true(1,numel(s));
idxs(dim) = false;
for i=1:numel(data)
  s_ = size(data{i});
  assert(all(s(idxs) == s_(idxs)));
end

if isempty(args.len)
  len = 1;
  for i=1:numel(data)
    len = lcm(len, size(data{i},dim));
  end
else
  len = args.len;
end
s(dim) = len;

out = data;
for i=1:numel(data)
  expand_factors = num2cell(ones(1,numel(s)));
  expand_factors{dim} = len / size(data{i},dim);
  out{i} = repelem(data{i}, expand_factors{:});
  out{i} = squeeze(out{i});
end
end
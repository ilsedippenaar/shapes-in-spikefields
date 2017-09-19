function out = normalizeData(data, varargin)
if isempty(data)
  out = [];
  return
end
if ~iscell(data)
  assert(isnumeric(data));
  data = {data};
else
  assert(numel(data) > 0);
end
p = inputParser;
p.addParameter('dim', find(size(data{1}) ~= 1, 1));
p.addParameter('len', []);
p.addParameter('method', []);
p.parse(varargin{:});
args = p.Results;

dim = args.dim;
s = size(data{1});
idxs = true(1,numel(s));
idxs(dim) = false;
for i=1:numel(data)
  s_ = size(data{i});
  assert(isequal(s(idxs), s_(idxs)), 'Dimensions along non-normalizing dimension do not match');
end

if isempty(args.len)
  len = 1;
  for i=1:numel(data)
    len = lcm(len, size(data{i},dim));
  end
else
  len = args.len;
end

out = data;
for i=1:numel(data)
  if isempty(args.method)
    expand_factors = num2cell(ones(1,ndims(data{i})));
    expand_factors{dim} = len / size(data{i},dim);
    out{i} = repelem(data{i}, expand_factors{:});
  else
    out{i} = apply(@(x) interp1(linspace(1,len,numel(x)),x,1:len,args.method),data{i},dim,len);
  end
  out{i} = squeeze(out{i});
end
end
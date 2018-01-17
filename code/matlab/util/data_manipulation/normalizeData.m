function out = normalizeData(data, varargin)
%NORMALIZEDATA Expands time series data to be the same length
% The data provided can be a numeric matrix or a cell array or matrices.
% This function accepts 3 keyword arguments. The dim (default first
% non-singleton dimension of the first element in data) sets the dimension
% to expand. len sets a size to expand to, and if unspecified will default
% to the least common multiple (LCM) of the lengths of the given data. This 
% can be very, very memory intensive if the lengths of the data given are 
% large and differ only slightly. The method parameter can be 'interp', in 
% which case linear interpolation is used to expand (useful with a set len
% parameter) or can be unspecified, in which case an LCM approach is taken.
%   For example: normalizeData({1:2,1:3}) returns {[1,1,1,2,2,2],
%   [1,1,2,2,3,3]}
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
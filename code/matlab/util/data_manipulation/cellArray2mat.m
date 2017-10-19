function out = cellArray2mat(c, data_type, fill)
if numel(c) == 0
  out = [];
  return
elseif numel(c) == 1
  out = c{1};
  return
end
assert(any(size(c) == 1));
dim = find(size(c) ~= 1);
assert(numel(dim) == 1 && (dim == 1 || dim == 2));
if nargin < 3
  if nargin < 2
    data_type = 'double';
  end
  if any(strcmp(data_type, {'double', 'single'}))
    fill = NaN;
  else
    fill = 0;
  end
end
if ~any(strcmp(data_type, {'double', 'single'})) && isnan(fill)
  error('Using data_type = %s with fill = NaN is not supported', data_type);
end

convertFunc = str2func(data_type);
if dim == 1
  out = repmat(convertFunc(fill), ...
    sum(cellfun(@(c) size(c,1), c)), max(cellfun(@(c) size(c,2),c)));
else
  out = repmat(convertFunc(fill), ...
    max(cellfun(@(c) size(c,1),c)), sum(cellfun(@(c) size(c,2), c)));
end
curr_idx = 1;
for i=1:numel(c)
  n = size(c{i},1);
  m = size(c{i},2);
  if dim == 1
    out(curr_idx:curr_idx+n-1,1:m) = convertFunc(c{i});
    curr_idx = curr_idx + n;
  else
    out(1:n,curr_idx:curr_idx+m-1) = convertFunc(c{i});
    curr_idx = curr_idx + m;
  end
end
end
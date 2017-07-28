function out = cellArray2mat(c, data_type)
if numel(c) == 0
  out = [];
  return
end
if nargin < 2
  data_type = 'double';
end
convertFunc = str2func(data_type);
if any(strcmp(data_type, {'double','single'}))
  out = nan(max(cellfun(@(c) size(c,1),c)), sum(cellfun(@(c) size(c,2), c)), data_type);
else
  out = zeros(max(cellfun(@(c) size(c,1),c)), sum(cellfun(@(c) size(c,2), c)), data_type);
end
curr_idx = 1;
for i=1:numel(c)
  n = size(c{i},1);
  m = size(c{i},2);
  out(1:n,curr_idx:curr_idx+m-1) = convertFunc(c{i});
  curr_idx = curr_idx + m;
end
end
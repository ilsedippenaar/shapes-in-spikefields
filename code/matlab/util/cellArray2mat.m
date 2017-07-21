function out = cellArray2mat(c, data_type)
if numel(c) == 0
  out = [];
  return
end
if nargin < 2
  data_type = 'double';
end
convertFunc = str2func(data_type);
out = zeros(size(c{1},1), sum(cellfun(@(c) size(c,2), c)), data_type);
curr_idx = 0;
for i=1:numel(c)
  n = size(c{1},2);
  out(:,curr_idx+1:curr_idx+n) = convertFunc(c{i});
  curr_idx = curr_idx + n;
end
end
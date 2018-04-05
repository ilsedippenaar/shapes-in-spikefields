function out = zipCell(varargin)
% makes the first level of c nested in the second level (like python's zip
% function)
if isempty(varargin)
  out = [];
  return
end
s = size(varargin{1});
assert(all(cellfun(@(d) all(size(d)==s), varargin)));
out = cell(s);
for i=1:numel(out)
  out{i} = cellfun(@(d) d{i}, varargin, 'UniformOutput', false);
end
end
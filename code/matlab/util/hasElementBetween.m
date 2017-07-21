function out = hasElementBetween(v, start, stop)
% tests whether v has an element such that it is in [start stop)
assert(stop >= start);
n = numel(v);
if n == 0
  out = false;
elseif v(1) >= start && v(end) < stop
  out = true;
else
  start_idx = binarySearch(v, start, ']');
  stop_idx = binarySearch(v, stop, '(');
  out = ~(isempty(start_idx) && isempty(stop_idx)) && (xor(isempty(start_idx), isempty(stop_idx)) || stop_idx - start_idx >= 0);
end
end
  
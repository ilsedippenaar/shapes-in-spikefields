function idx = binarySearch(v, x, mode, get_edges)
% Returns the idx such that x is in [v(idx) v(idx+1)).
% mode = '[' | '(' | ']' | ')'
% get_edges = logical
n = numel(v);
if nargin < 4
  get_edges = false;
end
if isempty(v) || (~get_edges  && (v(1) > x || v(end) < x))
  idx = [];
  return
elseif get_edges
  if v(1) > x
    if any(mode == '])')
      idx = 1;
    else
      idx = [];
    end
    return
  elseif v(end) < x
    if any(mode == '[(')
      idx = n;
    else
      idx = [];
    end
    return
  end
end

% short-circuiting is *very* important here
if any(mode == '])')
  condition = @(i) v(i) >= x && v(i) ~= x && v(i-1) < x || ...
                                v(i) == x && (i == n || v(i+1) > x);
else
  condition = @(i) v(i) <= x && v(i) ~= x && v(i+1) > x || ...
                                v(i) == x && (i == 1 || v(i-1) < x);
end

a = 1;
b = n;
idx = floor((a+b)/2);
while ~condition(idx)
  if v(idx) > x
    b = idx - 1;
  elseif v(idx) < x
    a = idx + 1;
  else
    if any(mode == '[(')
      b = idx - 1;
    else
      a = idx + 1;
    end
  end
  idx = floor((a+b)/2);
end

if v(idx) == x && mode == '('
  if idx == 1
    idx = [];
  else
    idx = idx - 1;
  end
  return
end
if v(idx) == x && mode == ')'
  if idx == n
    idx = [];
  else
    idx = idx + 1;
  end
  return
end
end

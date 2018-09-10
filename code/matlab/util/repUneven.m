function out = repUneven(v, n)
out = zeros(1,sum(n));
idx = 1;
for i=1:numel(n)
  out(idx:idx+n(i)-1) = v(i);
  idx = idx + n(i);
end
end
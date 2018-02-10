function union = unionInterval(a, b)
assert(a(2) >= a(1) && b(2) >= b(1));
if a(1) > b(2) || a(2) < b(1)
  union = [reshape(a,2,1), reshape(b,2,1)];
else
  union = [min(a(1),b(1)), max(a(2), b(2))];
end
end
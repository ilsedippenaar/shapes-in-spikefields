function test_binarySearch()
a = 1:10;
for i=a
  assert(binarySearch(a, i, '[') == i);
end
for i=a(1:end-1)
  assert(binarySearch(a, i, ')') == i+1);
end
assert(isempty(binarySearch(a, 10, ')')));

assert(binarySearch(a, 2.5, '[') == 2);
assert(binarySearch(a, 2.5, ')') == 3);

assert(binarySearch(a, 9.5, '[') == 9);
assert(binarySearch(a, 9.5, ')') == 10);

% outside
assert(isempty(binarySearch(a, 0, '[')));
assert(isempty(binarySearch(a, 0, ')')));

assert(isempty(binarySearch(a, 11, '[')));
assert(isempty(binarySearch(a, 11, ')')));

% odd-numbed array
b = [a 11];
for i=b
  assert(binarySearch(b, i, '[') == i);
end
for i=b(1:end-1)
  assert(binarySearch(b, i, ')') == i+1);
end
assert(isempty(binarySearch(b, 11, ')')));

% repetitions
c = [1 1 3 7 7 7 8 8 8];
assert(binarySearch(c, 7, '[') == 4);
assert(binarySearch(c, 7, ')') == 7);

assert(binarySearch(c, 1, '[') == 1);
assert(binarySearch(c, 1, ')') == 3);

assert(binarySearch(c, 8, '[') == 7);
assert(isempty(binarySearch(c, 8, ')'))); 


% silliness input
d = [];
assert(isempty(binarySearch(d, 1, '[')));
assert(isempty(binarySearch(d, 1, ')')));

e = 1;
assert(binarySearch(e, 1, '[') == 1);
assert(isempty(binarySearch(e, 1, ')')));

assert(isempty(binarySearch(e, 0, '[')));
assert(isempty(binarySearch(e, 0, ')')));

assert(isempty(binarySearch(e, 2, '[')));
assert(isempty(binarySearch(e, 2, ')')));

f = [1 2];
assert(binarySearch(f, 1, '[') == 1);
assert(binarySearch(f, 1, ')') == 2);

assert(binarySearch(f, 2, '[') == 2);
assert(isempty(binarySearch(f, 2, ')')));

fprintf('All tests passed\n');
end
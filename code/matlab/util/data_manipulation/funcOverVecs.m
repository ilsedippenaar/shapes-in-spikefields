function out = funcOverVecs(f, varargin)
assert(nargin(f) < 0 || nargin(f) == numel(varargin));
s = cellfun(@numel, varargin);
if numel(s) == 1
  s = [s,1];
end
out = zeros(s);
args = cartesianProd(varargin{:});
for i=1:prod(s)
  out(i) = f(args{i}{:});
end
end
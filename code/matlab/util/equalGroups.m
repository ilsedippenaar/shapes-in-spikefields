function varargout = equalGroups(v)
out = accumarray(v(:), 1:numel(v), [], @(x) {x});
m = min(cellfun(@numel, out));
varargout = cellfun(@(c) sort(datasample(c, m, 'Replace', false)), out, 'UniformOutput', false);
end
function varargout = funcOverCellArray(c, f, alloc_size, args, iter_args)
s = size(c);
varargout = cell(1,numel(alloc_size));
for i=1:numel(alloc_size)
  varargout{i} = zeros([alloc_size(i), s]);
end
for i=1:numel(c)
  out = f(c{i}, args{:}, iter_args{i}{:});
  for j=1:numel(varargout)
    varargout{j}(:,i) = out{i};
  end
end
for i=1:numel(varargout)
  varargout{i} = squeeze(varargout{i});
end
end
function out = funcInCells(data, func, empty_func, func_args, empty_func_args)
if nargin < 5 || isempty(empty_func_args)
  empty_func_args = {};
  if nargin < 4 || isempty(func_args)
    func_args = {};
    if nargin < 3 || isempty(empty_func)
      empty_func = @(x) [];
    end
  end
end
if iscell(data)
  out = data;
  for i=1:numel(data)
    out{i} = funcInCells(data{i}, func, empty_func, func_args, empty_func_args);
  end
elseif isempty(data)
  out = empty_func(empty_func_args{:});
else
  out = func(data, func_args{:});
end
end
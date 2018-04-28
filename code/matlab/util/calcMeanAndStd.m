function [m,s] = calcMeanAndStd(mat, dim, ignore_inf)
% robust against NaNs and any dimension
if nargin < 2
  dim = 1;
end
if nargin < 3
  ignore_inf = false;
end

if ignore_inf
  m = apply(@calcMeanInf, mat, dim, 1);
  s = apply(@calcStdInf, mat, dim, 1);
else
  m = nanmean(mat,dim);
  s = 1.96 * nanstd(mat,0,dim) ./ sqrt(sum(~isnan(mat),dim));
end

  function out = calcMeanInf(x)
    x = x(~isinf(x));
    x = x(:); % make a column
    out = calcMeanAndStd(x);
  end

  function out = calcStdInf(x)
    x = x(~isinf(x));
    x = x(:); % make column
    [~,out] = calcMeanAndStd(x);
  end
end
function plotMeanAndStds(mat, varargin)
if size(mat,2) > 0
  p = inputParser;
  p.addParameter('x', 1:size(mat,1));
  p.parse(varargin{:});
  args = p.Results;
  
  n = size(mat,2);
  if n == 1
    plot(args.x, mat);
    return
  end
  means = mean(mat, 2);
  stds = std(double(mat), 0, 2);
  plot(args.x, means, 'LineWidth', 2);
  hold on;
  plot(args.x, means + 1.96 * stds / sqrt(n), 'r:');
  plot(args.x, means - 1.96 * stds / sqrt(n), 'r:');
end
end
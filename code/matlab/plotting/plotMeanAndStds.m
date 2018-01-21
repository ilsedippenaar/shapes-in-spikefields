function [plt,x,y,std_err] = plotMeanAndStds(mat, varargin)
if size(mat,2) > 0
  p = inputParser;
  p.addParameter('x', 1:size(mat,1));
  p.parse(varargin{:});
  args = p.Results;
  
  x = args.x;
  n = size(mat,2);
  if n == 1
    plot(args.x, mat);
    x = args.x;
    y = mat;
    std_err = [];
    return
  end
  y = mean(mat, 2);
  stds = std(double(mat), 0, 2);
  plt = figure('Visible', 'off');
  plot(args.x, y, 'LineWidth', 2);
  hold on;
  std_err = 1.96 * stds / sqrt(n);
  plot(args.x, y + std_err, 'r:');
  plot(args.x, y - std_err, 'r:');
end
end
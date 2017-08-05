function out = padData(mat, to, data_type)
if nargin < 3
  data_type = 'double';
end
out_size = size(mat);
out_size(1) = to;
out = zeros(out_size, data_type);
out(1:size(mat,1), :) = mat;
end
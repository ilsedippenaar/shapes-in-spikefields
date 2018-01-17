function out = padData(mat, to, data_type)
%PADDATA Simple function to expand the first dimension of a matrix to a set
%length.
% The parameter mat should be a matrix, to specifies the length to pad to,
% and data_type (optional) sets the data type of the output matrix.
if nargin < 3
  data_type = 'double';
end
out_size = size(mat);
out_size(1) = to;
out = zeros(out_size, data_type);
out(1:size(mat,1), :) = mat;
end
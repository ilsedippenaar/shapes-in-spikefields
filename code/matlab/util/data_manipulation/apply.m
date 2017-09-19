function out = apply(f, data, dim, out_length)
if nargin < 4
  out_length = size(data,dim);
  if nargin < 3
    dim = 1;
  end
end

permute_vec = circshift(1:ndims(data), dim-1);
perm_data = permute(data, permute_vec);
perm_s = size(perm_data);
perm_out = zeros([out_length, perm_s(2:end)]);

data_idx = 1;
out_idx = 1;
for i=1:prod(perm_s(2:end))
  perm_out(out_idx:out_idx+out_length-1) = f(perm_data(data_idx:data_idx+perm_s(1)-1));
  data_idx = data_idx + perm_s(1);
  out_idx = out_idx + out_length;
end
out = ipermute(perm_out, permute_vec);
end
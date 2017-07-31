function mat = calculateAdjMatrix(vecs, dist_fun)
vecs = single(vecs);
if iscell(vecs)
  vecs = cellArray2mat(vecs, class(vecs{1}));
end
if nargin == 1
  dist_fun = @(v1,v2) (v1'*v2)/(norm(v1)*norm(v2));
elseif ischar(dist_fun)
  switch lower(dist_fun)
    case 'cosine'
      dist_fun = @(v1,v2) (v1'*v2)/(norm(v1)*norm(v2));
    case 'euclidean'
      dist_fun = @(v1,v2) norm(v1-v2);
    otherwise
      error('Unsupported function %s', dist_fun);
  end
end

n = size(vecs,2);
mat = zeros(n);
for i=1:n
  for j=i:n
    mat(i,j) = dist_fun(vecs(:,i),vecs(:,j));
    mat(j,i) = mat(i,j);
  end
end
end
function betweens = getBetweens(lfps, num_zeros)
if nargin < 2
  num_zeros = 3;
end
lfp = mean(lfps,2);

invalid_idxs = find(conv(lfp==0, ones(num_zeros,1), 'same') >= num_zeros)'; % 3 consecutive 0's = island of zeros
zero_betweens = getValidIntervals(invalid_idxs, num_zeros, numel(lfp));

lfp_limit = 2000;
high_window_size = 1000;
invalid_idxs = find(conv(abs(lfp) > lfp_limit, ones(high_window_size,1), 'same') >= 1)';
high_betweens = getValidIntervals(invalid_idxs, high_window_size, numel(lfp));

% intersect zero_betweens and high_betweens
betweens = cell(1,size(zero_betweens,2)*size(high_betweens,2));
for i=1:size(zero_betweens,2)
  for j=1:size(high_betweens,2)
    betweens{i*(size(high_betweens,2)-1)+j} = intersectInterval(zero_betweens(:,i), high_betweens(:,j))';
  end
end
betweens = [betweens{:}];

  function valid_intervals = getValidIntervals(invalid_idxs, window_size, total_size)
    if isempty(invalid_idxs)
      valid_intervals = [1; total_size];
      return
    end
    edges = find(diff(invalid_idxs)>1);
    beg_offset = ceil((window_size+1)/2); % convolution silliness and indexing - odd length window is symmetric, even isn't
    end_offset = ceil(window_size/2);
    valid_sections = [1, invalid_idxs(edges)+beg_offset, invalid_idxs(end)+beg_offset;
                     invalid_idxs(1)-end_offset, invalid_idxs(edges+1)-end_offset, total_size];
    valid_intervals = valid_sections(:, diff(valid_sections,1) > 0); % check for any invalid sections at the beginning or end
  end
end
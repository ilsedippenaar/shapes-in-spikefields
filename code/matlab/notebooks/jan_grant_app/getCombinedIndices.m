function indices = getCombinedIndices(lfps)
% lfps is cell array with days -> electrodes -> time x trials int16 array

% make sure everything is the same size
num_electrodes = numel(lfps{1});
for i=1:numel(lfps)
  assert(num_electrodes == numel(lfps{i}));
end

% number of trials in a day = maximum number of trials for an electrode
total_in_days = cell(numel(lfps),1);
for i=1:numel(lfps)
  num_in_day = cell(1,num_electrodes);
  for j=1:num_electrodes
    num_in_day{j} = ones(size(lfps{i}{j},2),1,'logical');
  end
  % trials in day x num_electrodes
  total_in_days{i} = cellArray2mat(num_in_day, 'logical', 0);
end
% total trials x num_electrodes
total_in_days = cellArray2mat(total_in_days);

indices = cell(1,num_electrodes);
for i=1:num_electrodes
  indices{i} = find(total_in_days(:,i));
end
end
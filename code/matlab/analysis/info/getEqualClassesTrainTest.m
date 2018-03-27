function [train,test] = getEqualClassesTrainTest(d)
labels = categorical(d.label);
n = min(countcats(labels));
labs = unique(labels);
num_labs = numel(labs);
assert(num_labs == 2);
train = cell(num_labs,1);
test = cell(num_labs,1);
for i=1:num_labs
  lab = labs(i);
  subset = d(labels == lab,:);
  idxs = randperm(size(subset, 1));
  shuffled = subset(idxs,:);
  % get 50% of both categories so that fitting function isn't biased
  % towwards one class
  train{i} = shuffled(1:ceil(n/2), :); 
  test{i} = shuffled(end-ceil(n/2):end, :);
end
train = vertcat(train{:});
test = vertcat(test{:});
end
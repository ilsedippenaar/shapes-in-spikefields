function labs = getLabels(tbl, name, num_splits)
data = tbl.(name);
bins = quantile(data, (0:num_splits)/num_splits, 1);
labs = apply(@(i) discretize(data(:,i), bins(:,i)), ...
                 1:size(data,2), 1, size(data,1)); % discretize each row (row = trial)
end
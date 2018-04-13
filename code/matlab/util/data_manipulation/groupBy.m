function [groups,vals] = groupBy(tbl, var)
[G, vals] = findgroups(tbl.(var));
groups = cell(numel(vals),1);
for i=1:numel(groups)
  groups{i} = tbl(i==G, :);
end
end
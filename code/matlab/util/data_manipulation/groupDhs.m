function grouped_dhs = groupDhs(dhs)
days = categorical({dhs.date});
unique_days = categories(days);
grouped_dhs = cell(1,numel(unique_days));
for i=1:numel(unique_days)
  grouped_dhs{i} = dhs(unique_days(i) == days);
end
end
function out = clean_lfp_selection(lfp_selection)
first_not_empty = find(cellfun(@numel, lfp_selection{1}),1);
last_not_empty = find(cellfun(@numel, lfp_selection{1}),1, 'last');
out = lfp_selection;
for i=1:numel(lfp_selection)
  out{i} = lfp_selection{i}(first_not_empty:last_not_empty);
end
end
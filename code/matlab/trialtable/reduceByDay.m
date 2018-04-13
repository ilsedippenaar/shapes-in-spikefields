function reduceByDay(intermKey, intermValIter, outKVStore)
if strcmp(intermKey, 'coherence')
  out = {};
  while hasnext(intermValIter)
    % average across trials
    cohs = nanmean(getnext(intermValIter),4);
    out{end+1} = cohs;
  end
  out = cat(4, out{:});
  add(outKVStore, intermKey, out);
end
end
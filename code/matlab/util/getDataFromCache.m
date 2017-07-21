function data = getDataFromCache(cache_dir, params)
if exist(fullfile(cache_dir, 'cache_map.mat'), 'file') ~= 2
  data = [];
else
  cache_map_file = load(fullfile(cache_dir, 'cache_map.mat'));
  cache_map = cache_map_file.cache_map;
  data = [];
  for i=1:numel(cache_map)
    if isequal(cache_map(i).params, params)
      data_file = load(fullfile(cache_dir, cache_map(i).filename));
      data = data_file.data;
      break
    end
  end
end
end
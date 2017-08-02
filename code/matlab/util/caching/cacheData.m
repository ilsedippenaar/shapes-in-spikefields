function cacheData(cache_dir, params, data)
if exist(fullfile(cache_dir, 'cache_map.mat'), 'file') ~= 2
  cache_map = [];
else
  cache_map_file = load(fullfile(cache_dir, 'cache_map.mat'));
  cache_map = cache_map_file.cache_map;
end
n = numel(cache_map);
cache_map(n+1).params = params;
cache_map(n+1).filename = [num2str(n+1) '.mat'];
save(fullfile(cache_dir, 'cache_map.mat'), 'cache_map');
save(fullfile(cache_dir, cache_map(end).filename), 'data');
end
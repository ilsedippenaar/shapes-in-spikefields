function [cohs, freqs] = calculateCoherence(lfp_slices, cache_dir, conditions, varargin)
p = inputParser;
p.addParameter('Fs', 1000);
p.addParameter('window_width', 512);
p.addParameter('num_fft', 512);
p.parse(varargin{:});
args = p.Results;

if ~mod(args.num_fft, 2)
  num_points = args.num_fft/2 + 1;
else
  num_points = (args.num_fft+1)/2;
end
freqs = [];
n = numel(lfp_slices);
cohs = zeros(num_points, n, n);
params = args;
params.conditions = conditions;
params.type = 'coherence';
cached_data = getDataFromCache(cache_dir, params);

if isempty(cached_data)
  for i=1:n
    fprintf('Electrode: %d\n', i);
    r = int32(rem(numel(lfp_slices{i}), args.window_width)); % this is all to prevent what is presumably a bug in mscohere
    cohs(:,i,i) = 1;
    for j=i+1:n
      [cxy,f] = mscohere(lfp_slices{i}(1:end-r,:), lfp_slices{j}(1:end-r, :), hann(args.window_width),...
        args.window_width/2, args.num_fft, args.Fs);
      assert(isempty(freqs) || all(freqs == f));
      freqs = f;
      cohs(:,i,j) = mean(cxy,2);
      cohs(:,j,i) = cohs(:,i,j);
    end
  end
  cached_data = [];
  cached_data.cohs = cohs;
  cached_data.freqs = freqs;
  cacheData(cache_dir, params, cached_data);
else
  fprintf('Found cached coherences\n');
  cohs = cached_data.cohs;
  freqs = cached_data.freqs;
end
end
function [psds, freqs, pconf] = calculatePsd(data, cache_dir, conditions, varargin)
p = inputParser;
p.addParameter('method', 'welch');
p.addParameter('Fs', 1000);
p.addParameter('window_width', 512); % for welch
p.addParameter('time_halfbandwidth', 2.5); % for mtm
p.addParameter('num_fft', 512);
p.parse(varargin{:});
args = p.Results;

if ~mod(args.num_fft, 2)
  num_points = args.num_fft/2 + 1;
else
  num_points = (args.num_fft+1)/2;
end
freqs = [];
n = numel(data);
psds = zeros(num_points, n);

% get from cache if available
params = args;
if strcmp(args.method, 'pwelch')
  params.time_halfbandwidth = [];
else
  params.window_width = [];
end
params.conditions = conditions;
params.type = 'psd';
cached_data = getDataFromCache(cache_dir, params);
if isempty(cached_data)
  for i=1:n
%     r = int32(rem(numel(data{i}), args.window_width));
    if strcmp(args.method, 'pwelch')
      [psd,f,pconf] = pwelch(data{i}, hann(args.window_width),...
        args.window_width/2, args.num_fft, args.Fs, 'ConfidenceLevel', 0.95);
    else
      [psd,f,pconf] = pmtm(data{i}, args.time_halfbandwidth, args.num_fft, args.Fs, 'ConfidenceLevel', 0.95);
    end
    assert(isempty(freqs) || all(freqs == f));
    freqs = f;
    psds(:,i) = mean(psd,2);
  end
  cached_data = [];
  cached_data.psds = psds;
  cached_data.freqs = freqs;
  cached_data.pconf = pconf;
  cacheData(cache_dir, params, cached_data);
else
  fprintf('Found cached PSDs\n');
  psds = cached_data.psds;
  freqs = cached_data.freqs;
  pconf = cached_data.pconf;
end
end
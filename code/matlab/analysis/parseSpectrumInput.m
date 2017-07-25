function [cached_data, cache_params, out, args, is_variable_len, varargout] = parseSpectrumInput(type, in_data, cache_dir, conditions, varargin)
p = inputParser;
p.addParameter('method', 'welch');
p.addParameter('Fs', 1000);
p.addParameter('window_width', 512);
p.addParameter('num_fft', 512); % only NOT used if method ~= welch and is_variable_len is true
p.addParameter('T', 0.5);
p.addParameter('W', 6);
p.parse(varargin{:});
args = p.Results;
is_variable_len = iscell(in_data{1});

% get data from cache (and create dpss)
cache_params = [];
cache_params.type = type;
cache_params.method = args.method;
cache_params.Fs = args.Fs;
cache_params.is_variable_len = is_variable_len;
cache_params.conditions = conditions;
if strcmp(args.method, 'welch')
  cache_params.window_width = args.window_width;
  cache_params.num_fft = args.num_fft;
else
  cache_params.T = args.T;
  cache_params.W = args.W;
end
cached_data = getDataFromCache(cache_dir, cache_params);

% allocate output
n = numel(in_data);
if strcmp(args.method, 'welch')
  if ~mod(args.num_fft, 2)
    num_points = args.num_fft/2 + 1;
  else
    num_points = (args.num_fft+1)/2;
  end
else
  mt_params = [];
  mt_params.Fs = args.Fs;
  mt_params.pad = 1;
  mt_params.trialave = true;
  % precalculate dpss tapers
  % TODO: use native matlab implementation of dpss 
  if is_variable_len
    varargout{2} = arrayfun(@(x) dpsschk([args.T*args.W, 2*args.T*args.W-1], x, args.Fs), ...
      cellfun(@numel, in_data{1}), 'UniformOutput', false);
    num_points = 2^(nextpow2(args.num_fft)-1)+1;
  else
    mt_params.tapers = dpsschk([args.T*args.W, 2*args.T*args.W-1], args.T, args.Fs);
    num_points = 2^nextpow2(size(in_data{1},1)-1)+1;
  end
  varargout{1} = mt_params;
end
if strcmp(type, 'psd')
  out = zeros(num_points, n);
else
  out = zeros(num_points, n, n);
end
end
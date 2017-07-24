function [cohs, freqs] = calculateCoherence(lfp_slices, cache_dir, conditions, varargin)
p = inputParser;
p.addParameter('method', 'welch');
p.addParameter('Fs', 1000);
p.addParameter('window_width', 512);
p.addParameter('num_fft', 512); % only NOT used if method ~= welch and lfp_slices is a cell array of matrices
p.addParameter('T', 0.5);
p.addParameter('W', 6);
p.parse(varargin{:});
args = p.Results;
is_variable_len = iscell(lfp_slices{1});

% get data from cache (and create dpss)
cache_params = [];
cache_params.type = 'coherence';
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
  mt_params = [];
  mt_params.Fs = args.Fs;
  mt_params.pad = 1;
  mt_params.trialave = true;
  % precalculate dpss tapers
  % TODO: use native matlab implementation of dpss 
  if is_variable_len
    tapers = arrayfun(@(x) dpsschk([args.T*args.W, 2*args.T*args.W-1], x, args.Fs), ...
      cellfun(@numel, lfp_slices{1}), 'UniformOutput', false);
  else
    mt_params.tapers = dpsschk([args.T*args.W, 2*args.T*args.W-1], args.T, args.Fs);
  end
end
cached_data = getDataFromCache(cache_dir, cache_params);

% allocate output
n = numel(lfp_slices);
if strcmp(args.method, 'welch')
  if ~mod(args.num_fft, 2)
    num_points = args.num_fft/2 + 1;
  else
    num_points = (args.num_fft+1)/2;
  end
  cohs = zeros(num_points, n, n);
else
  if is_variable_len
    cohs = zeros(2^nextpow2(args.num_fft)+1,n,n);
  else
    cohs = zeros(2^nextpow2(size(lfp_slices{1},1))+1,n,n);
  end
end
freqs = [];

% TODO: refactor all of this
% calculate coherence
if isempty(cached_data)
  for i=1:n
    fprintf('Electrode: %d\n', i);
    cohs(:,i,i) = 1;
    for j=i+1:n
      if is_variable_len
        m = 0;
        for k=1:numel(lfp_slices{i})
          if strcmp(args.method, 'welch')
            [cxy,f] = mscohere(lfp_slices{i}{k}, lfp_slices{j}{k}, ...
              hann(args.window_width), args.window_width/2, args.num_fft, args.Fs);
            cohs(:,i,j) = cohs(:,i,j) + cxy;
            m = m + 1;
          else
            if 2^nextpow2(lfp_slices{i}{k}) < args.num_fft
              mt_params.tapers = tapers{k};
              [cxy,~,~,~,~,f] = coherencyc(lfp_slices{i}, lfp_slices{j}, mt_params);
              factor = numel(cxy)-1 / arg.num_fft;
              if factor ~= 1
                cxy = downsample(cxy, factor);
                f = downsample(f, factor);
              end
              cohs(:,i,j) = cxy;
              m = m + 1;
            end
          end
        end
        cohs(:,i,j) = cohs(:,i,j) / m;
      else
        %r = int32(rem(size(lfp_slices{i},1), args.window_width)); % this is all to prevent what is presumably a bug in mscohere
        if strcmp(args.method, 'welch')
          [cxy,f] = mscohere(lfp_slices{i}, lfp_slices{j}, ...
            hann(args.window_width), args.window_width/2, args.num_fft, args.Fs);
          cohs(:,i,j) = mean(cxy,2);
        else
          % TODO use error outputs
          [cxy,~,~,~,~,f] = coherencyc(lfp_slices{i}, lfp_slices{j}, mt_params);
          cohs(:,i,j) = cxy;
        end
      end
      assert(isempty(freqs) || all(freqs == f));
      freqs = f;
      cohs(:,j,i) = cohs(:,i,j);
    end
  end
  cached_data = [];
  cached_data.cohs = cohs;
  cached_data.freqs = freqs;
  cacheData(cache_dir, cache_params, cached_data);
else
  fprintf('Found cached coherences\n');
  cohs = cached_data.cohs;
  freqs = cached_data.freqs;
end
end
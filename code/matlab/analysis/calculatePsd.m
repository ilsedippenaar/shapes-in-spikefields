function [psds, freqs] = calculatePsd(lfps, cache_dir, conditions, varargin)
[cached_data, cache_params, psds, args, is_variable_len, mt_params, tapers] = parseSpectrumInput('psd', lfps, cache_dir, conditions, varargin{:});
freqs = [];

if isempty(cached_data)
  n = size(psds,2);
  num_points = size(psds,1);
  for i=1:n
    if is_variable_len
      m = 0;
      for j=1:numel(lfps{i})
        if strcmp(args.method, 'welch')
          [psd,f] = pwelch(lfps{i}{j},...
            hann(args.window_width), args.window_width/2, args.num_fft, args.Fs);
          psds(:,i) = psds(:,i) + psd;
          m = m + 1;
        else
          % TODO: is this condition necessary? Shouldn't very small samples
          % be excluded elsewhere if at all? It would also make this
          % implementation much cleaner and generic
          if nextpow2(numel(lfps{i}{j})) >= nextpow2(args.num_fft)
            mt_params.tapers = tapers{j};
            [psd,f] = mtspectrumc(single(lfps{i}{j}), mt_params);
            factor = (numel(psd)-1) / (num_points-1);
            if factor ~= 1
              psd = downsample(psd, factor);
              f = downsample(f, factor);
            end
            psds(:,i) = psds(:,i) + psd;
            m = m + 1;
          end
        end
      end
      psds(:,i) = psds(:,i) / m;
    else
      if strcmp(args.method, 'welch')
        [psd,f] = pwelch(lfps{i}, ...
          hann(args.window_width), args.window_width/2, args.num_fft, args.Fs);
        psds(:,i) = mean(psd,2);
      else
        % TODO use error outputs
        [psd,f] = mtspectrumc(lfps{i}, mt_params);
        psds(:,i) = psd;
      end
    end
    assert(isempty(freqs) || all(freqs == f));
    freqs = f;
  end
  cached_data = [];
  cached_data.psds = psds;
  cached_data.freqs = freqs;
  cacheData(cache_dir, cache_params, cached_data);
else
  fprintf('Found cached psds\n');
  psds = cached_data.psds;
  freqs = cached_data.freqs;
end
end
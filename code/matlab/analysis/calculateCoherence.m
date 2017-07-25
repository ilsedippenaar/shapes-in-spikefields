function [cohs, freqs] = calculateCoherence(lfps, cache_dir, conditions, varargin)
[cached_data, cache_params, cohs, args, is_variable_len, mtm_objs] = parseSpectrumInput('coherence', lfps, cache_dir, conditions, varargin{:});
if strcmp(args.method, 'mtm')
  mt_params = mtm_objs{1};
  if is_variable_len
    tapers = mtm_objs{2};
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
        for k=1:numel(lfps{i})
          if strcmp(args.method, 'welch')
            [cxy,f] = mscohere(lfps{i}{k}, lfps{j}{k}, ...
              hann(args.window_width), args.window_width/2, args.num_fft, args.Fs);
            cohs(:,i,j) = cohs(:,i,j) + cxy;
            m = m + 1;
          else
            if nextpow2(numel(lfps{i}{k})) >= nextpow2(args.num_fft)
              mt_params.tapers = tapers{k};
              [cxy,~,~,~,~,f] = coherencyc(single(lfps{i}{k}), single(lfps{j}{k}), mt_params);
              factor = (numel(cxy)-1) / (num_points-1);
              if factor ~= 1
                cxy = downsample(cxy, factor);
                f = downsample(f, factor);
              end
              cohs(:,i,j) = cohs(:,i,j) + cxy;
              m = m + 1;
            end
          end
        end
        cohs(:,i,j) = cohs(:,i,j) / m;
      else
        %r = int32(rem(size(lfp_slices{i},1), args.window_width)); % this is all to prevent what is presumably a bug in mscohere
        if strcmp(args.method, 'welch')
          [cxy,f] = mscohere(lfps{i}, lfps{j}, ...
            hann(args.window_width), args.window_width/2, args.num_fft, args.Fs);
          cohs(:,i,j) = mean(cxy,2);
        else
          % TODO use error outputs
          [cxy,~,~,~,~,f] = coherencyc(lfps{i}, lfps{j}, mt_params);
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
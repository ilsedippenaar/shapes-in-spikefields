function [smoothed,t] = calculateSmoothedSpikes(spikes, window, start, len)
if isnumeric(spikes)
  if size(spikes,1) == 1
    spikes = spikes';
  end
  spikes = {spikes};
end
spikes = expandSpikes(spikes, start, len);
out_len = len - numel(window) + 1;
smoothed = apply(@(x) conv(double(x), window, 'valid'), spikes, 1, out_len);
offset = (numel(window)-1)/2;
t = start+offset:start+len-1-offset;
end
function [smoothed,t] = calculateSmoothedSpikes(spikes, window, start, len)
if isnumeric(spikes)
  if size(spikes,1) == 1
    spikes = spikes';
  end
  spikes = {spikes};
end
spikes = expandSpikes(spikes, start, len);
smoothed = apply(@(x) conv(double(x), window, 'same'), spikes, 1, len);
t = start:start+len-1;
end
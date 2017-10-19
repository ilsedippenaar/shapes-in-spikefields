function [plts,psth,t] = plotPsth(spikes, window_size, start, len)
if isnumeric(spikes)
  if size(spikes,1) == 1
    spikes = spikes';
  end
  spikes = {spikes};
end

plts = gobjects(1,numel(spikes));
for i=1:numel(spikes)
  plts(i) = figure('Visible', 'off');
  histogram(spikes{i},start:window_size:start+len);
  xlabel('Time (ms)');
  ylabel('Spike count');
end
psth = cellArray2mat(cellfun(@(c) histcounts(c, start:window_size:start+len)', spikes, 'UniformOutput', false));
t = start:window_size:start+len;
end
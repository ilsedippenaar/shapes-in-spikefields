function out = expandSpikes(spikes, start, len)
if ~iscell(spikes)
  spikes = {spikes};
end
out = zeros(len,numel(spikes),'int8');
for i=1:numel(spikes)
  if ~isempty(spikes{i})
    assert(max(spikes{i}) <= start+len-1);
    out(spikes{i}-start+1,i) = 1;
  end
end
end
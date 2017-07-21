function [coeffs,S,mu] = regressReactionTimes(dh, lfps, start_time)
mapping = [repmat({'fixate'},1,numel(dh.fixate)), ...
  repmat({'noise'},1,numel(dh.noise)), ...
  repmat({'shape'},1,numel(dh.shape)), ...
  repmat({'saccade'},1,numel(dh.saccade))];
[all_events, idxs] = sort([dh.fixate, dh.noise, dh.shape, dh.saccade]);
labels = mapping(idxs);
labels = labels(all_events > start_time);
all_events = all_events(all_events > start_time);
rxn_times = [];
slopes = [];
for i=2:numel(all_events)-1
  if strcmp(labels{i}, 'shape') && strcmp(labels{i+1}, 'saccade') && all_events(i)-all_events(i-1) > 500
    rxn_time = double(all_events(i+1)-all_events(i));
    y = mean(lfps(all_events(i)-200:all_events(i+1),:),2);
    coeffs = [ones(numel(y),1) (0:numel(y)-1)'] \ y; % taking slope of the mean (not mean of the slopes)
    rxn_times = [rxn_times rxn_time];
    slopes = [slopes coeffs(2)];
  end
end

[coeffs,S,mu] = polyfit(slopes, rxn_times, 1);
csvwrite('slope_data.csv', [slopes' rxn_times']);
end
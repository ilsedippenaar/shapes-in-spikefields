function plt = plotFullTraces(lfps, elec_num)
num_in_lfp = cellfun(@numel, lfps{elec_num});
n = max(num_in_lfp);
mean_noise = zeros(n,1);
std_noise = zeros(n,1);
for i=1:n
  vals = [];
  for j=1:numel(lfps{elec_num})
    if ~isempty(lfps{elec_num}{j}) && numel(lfps{elec_num}{j}) >= i
      vals = [vals lfps{elec_num}{j}(i)];
    end
  end
  if ~isempty(vals)
    mean_noise(i) = mean(vals);
    std_noise(i) = std(single(vals));
  end
end
first_not_empty = find(num_in_lfp, 1);
last_not_empty = find(num_in_lfp, 1, 'last');
num_range = zeros(1,n);
num_range(num_in_lfp(first_not_empty:last_not_empty)) = 1;
num_counted = cumsum(num_range, 2, 'reverse');

plt = figure('Visible', 'off');
subplot(2,1,1);
plot(num_counted(num_counted > 1));
title('Number of trial counted at given time')

subplot(2,1,2);
plot(mean_noise(num_counted > 1));
hold on
plot(mean_noise(num_counted > 1) + std_noise(num_counted > 1) * 1.96 *[1,-1], 'r--');
title(sprintf('LFP trace for electrode %d', elec_num));
end
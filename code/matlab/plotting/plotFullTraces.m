function [plt,mean_lfp,std_lfp] = plotFullTraces(lfps, elec_num)
[mean_lfp,std_lfp,num_counted] = combineVariableLengthLfps(lfps);

plt = figure('Visible', 'off');
subplot(2,1,1);
plot(num_counted(num_counted > 1));
title('Number of trial counted at given time')

subplot(2,1,2);
plot(mean_lfp(num_counted > 1));
hold on
plot(mean_lfp(num_counted > 1) + std_lfp(num_counted > 1) * 1.96 *[1,-1], 'r--');
title(sprintf('LFP trace for electrode %d', elec_num));
end
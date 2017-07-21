function test_getLfpSlices(dh, between)
selection = {'between', between, 'melt', true, 'type', 'lfp'};
a=dh.getLfpSlices('0,any,0,128,saccade,1,0', between(1), between(2), 'selection_params', selection);

saccade = int32([dh.trials.saccade]);
test_lfps = (1:saccade(5)+101)';
b = dh.getLfpSlices('0,saccade,1,100,any,0,0', between(1), between(2), 'lfps', test_lfps);
c = dh.getLfpSlices('0,any,0,1,saccade,1,100', between(1), between(2), 'lfps', test_lfps);
d = dh.getLfpSlices('0,saccade,10000,1,saccade,1,100', between(1), between(2), 'lfps', test_lfps);
e = dh.getLfpSlices('0,any,0,10000,any,0,0', between(1), between(2), 'lfps', test_lfps);
f = dh.getLfpSlices('100,saccade,100,0,any,0,0', between(1), between(2), 'lfps', test_lfps);

fprintf('All passed!\n');
end
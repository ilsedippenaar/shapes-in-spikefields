function test_getDataSlices(dh)
conditions(1).vec = [5 20 53 80 91];
conditions(1).negate = false;
conditions(1).range = [-1 0];
data = {1:100, 4:103};
slices = dh.getDataSlices(data, 'lfp', [-20 10], conditions);
assert(size(slices{1},2) == 2);
assert(all(slices{1}(20,:) == [20 53]));

slices = dh.getDataSlices(data, 'lfp', [-19 10], conditions);
assert(size(slices{1},2) == 2);
assert(all(slices{1}(19,:) == [20 53]));

slices = dh.getDataSlices(data, 'lfp', [-21 10], conditions);
assert(size(slices{1},2) == 1);
assert(all(slices{1}(21,:) == [53]));

conditions(2).vec = [17];
conditions(2).negate = true;
conditions(2).range = [-10 10];
slices = dh.getDataSlices(data, 'lfp', [-5 5], conditions);
assert(size(slices{1},2) == 2);
assert(all(slices{1}(5,:) == [53 80]));

fprintf('All passed\n');
end
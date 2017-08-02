if multiple_days && exist('dhs', 'var') ~= 1
  dhs = DataHandler.fromDates(datetime('2013-05-1'):datetime('2013-05-29'), ...
    fullfile(data_dir, 'trials'), ...
    fullfile(data_dir, 'lfps'), ...
    fullfile(data_dir, 'data_handlers'), 'clean', true);
end

betweens = {
  [
    1, 1.254e6;
    1.444e6, 4.193e6
  ]', ...
  [
    2.163e4, 1.092e6;
    1.233e6, 2.851e6;
    3.196e6, 3.889e6
  ]', ...
  [
    1, 1.265e6;
    1.302e6, 2.435e6;
    2.493e6, 3.091e6;
    3.129e6, 3.863e6
  ]', ...
  [
    1, 2.361e6;
    2.577e6, 3.779e6;
    3.803e+6, 4.088e+6
  ]', ...
  [
    1, 1.314e6;
    1.862e+6, 2.609e+6;
    2.639e+6, 3.322e+6;
    3.364e+6, 3.511e+6;
    3.538e+6, 5.161e+6
  ]', ...
  [
    1, 1.130e+6;
    1.177e+6, 1.457e+6;
    1.744e+6, 1.982e+6
  ]', ...
  [
    5.591e+4, 1.487e+6
  ]', ...
  [
    1.566e+4, 3.977e+5;
    4.081e+5, 9.102e+5;
    9.308e+5, 1.007e+6
  ]', ...
  [
    1.113e+5, 3.098e+6
  ]', ...
  [
    1, 3.803e+6;
    3.836e+6, 4.561e+6
  ]', ...
  [
    1, 4.869e+5;
    5.107e+5, 1.043e+6;
    1.067e+6, 3.713e+6
  ]', ...
  [
    4.134e+4, 2.112e+6;
    2.142e+6, 2.593e+6;
    2.675e+6, 3.641e+6
  ]', ...
  [
    2.159e+5, 2.370e+6;
    2.415e+6, 3.953e+6;
    4.224e+6, 4.776e+6
  ]', ...
  [
    1, 4.512e+6
  ]', ...
  [
    1, 2.545e+6;
    2.979e+6, 4.115e+6
  ]', ...
  [
    1, 1.778e+6;
    1.808e+6, 2.329e+6;
    2.364e+6, 3.255e+6;
    3.776e+6, 4.538e+6
  ]', ...
  [
    1, 3.348e+5;
    3.456e+5, 8.402e+5
  ]', ...
  [
    1, 1.226e+6;
    1.254e+6, 1.742e+6;
    1.767e+6, 2.258e+6;
    2.292e+6, 3.351e+6
  ]'};
all_dhs = cell(1,sum(cellfun(@(c) size(c,2), betweens)));
idx = 1;
for j=1:numel(betweens)
  num = size(betweens{j},2);
  all_dhs(idx:idx+num-1) = dhs(j).split(betweens{j});
  idx = idx + num;
end
for j=1:numel(all_dhs)
  if all_dhs{j}.num_trials == 0
    all_dhs{j} = [];
  end
end
dhs = [all_dhs{:}];
clear all_dhs
%%
trial_sections = [4,5];
stim_lfps = cell(96,numel(trial_sections));
for i=1:numel(trial_sections)
  all_stim = cell(1,numel(dhs));
  for j=1:numel(dhs)
    stim = dhs(j).select('type', 'lfp', 'trial_section', trial_sections(i), 'trial_result', {'true_positive'});
    for k=1:numel(stim)
      stim{k} = cellArray2mat(stim{k}','single');
    end
    all_stim{j} = expandCellArray(dhs(j), stim);
  end
  stim_lfps(:,i) = combineCellArrays('single', all_stim{:});
end
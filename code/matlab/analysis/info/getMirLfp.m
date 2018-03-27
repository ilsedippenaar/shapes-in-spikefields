function mir = getMirLfp(d)
% d is a table with columns 'data' and 'label'
if isempty(d)
  mir = [];
  return
end
[train,test] = getEqualClassesTrainTest(d);

% fit and evaluate the model
mdl = fitcdiscr(single(train.data), train.label);
pred_labs = categorical(predict(mdl, single(test.data)));
mi = mutInfo(double(pred_labs), double(categorical(test.label)));
mir = mi / size(d.data,2) * 1000; % MI per second
end
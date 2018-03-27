function [mir, se, percentile] = calculateInfo(dhs, data_type, info_type, binwidth, delay, varargin)
p = inputParser;
p.addParameter('length_postnoise_response', 250);
p.addParameter('average', true);
p.addParameter('expand_electrodes', true);
p.addParameter('single_unit', true);
p.addParameter('quantile', 0.95);
p.parse(varargin{:});
args = p.Results;

% validate input
assert(isa(dhs, 'DataHandler'));
if strcmp(data_type, 'spike')
  assert(numel(unique([dhs.num_units])) == 1);
elseif numel(dhs) > 1
  % check for different electrode mappings, if there are, then set expand_electrodes to true
  maps = [dhs.electrode_mapping];
  electrodes = arrayfun(@(x) maps(:,x), 1:2:size(maps,2), 'UniformOutput', false);
  if ~isequal(electrodes{:})
    args.expand_electrodes = true;
  end
end
if args.average % averaging takes priority over expand electrodes
  args.expand_electrodes = false;
end
if isempty(data_type)
  data_type = '';
end
validatestring(data_type, {'spike', 'lfp', ''});
validatestring(info_type, {'sensory', 'choice', 'behavioral'});

% define selection parameters
if strcmp(info_type, 'sensory')
  conds = {
    cell2struct([...
      {'noise', false, [0,1]}; ...
      {'shape', true, [0, delay+binwidth]}], {'vec','negate','range'}, 2), ...
    cell2struct([...
      {'shape', false, [0,1]}; ...
      {'noise', true, [-args.length_postnoise_response, delay+binwidth]}], {'vec','negate','range'}, 2)};
    labels = {'noise', 'shape'};
    select_ranges = {delay+[0,binwidth], delay+[0,binwidth]};
elseif strcmp(info_type, 'choice')
  conds = {
    cell2struct([...
      {'fixate', false, [0,1]}; ...
      {'saccade', true, [0, delay+binwidth]}], {'vec', 'negate', 'range'}, 2), ...
    cell2struct([...
      {'saccade', false, [0,1]}; ...
      {'fixate', true, [-delay-binwidth, 0]}], {'vec', 'negate', 'range'}, 2)};
    labels = {'fixate', 'saccade'};
    select_ranges = {delay+[0,binwidth], [-binwidth,0]-delay};
else
  mir = getBehavioralInfo(dhs, binwidth, delay);
  return
end

% select data
if strcmp(data_type, 'spike')
  n = dhs(1).num_units;
else
  if args.average
    n = 1;
  elseif args.expand_electrodes
    n = 96;
  else
    n = dhs(1).num_lfp_electrodes;
  end
end
data = cell(1,numel(conds));
for i=1:numel(conds)
  data{i} = cell(n, numel(dhs));
  for j=1:numel(dhs)
    d = dhs(j).getDataSlices(data_type, select_ranges{i}, conds{i});
    if strcmp(data_type, 'lfp') && args.average
      mat_lfp = reshape([d{:}], [size(d{1}), numel(d)]); % time x trials x electrodes
      data{i}{1,j} = mean(mat_lfp, 3);
    elseif strcmp(data_type, 'lfp') && args.expand_electrodes
      data{i}(:,j) = expandCellArray(dhs(j).electrode_mapping, d);
    else
      data{i}(:,j) = d;
    end
  end
end
% check that there are equal number of units / electrodes / averaged electrodes
assert(numel(unique(cellfun(@(c) size(c,1), data))) == 1);
n = size(data{1},1);

% build tables for analysis and calculate MI rate
mir = cell(1,n);
se = cell(size(mir));
percentile = cell(size(mir));
for i=1:n
  combined = cell(1,numel(data));
  for j=1:numel(data) % conds
    d = [data{j}{i,:}]';
    combined{j} = table(d, repelem(labels(j), size(d,1), 1), 'VariableNames', {'data', 'label'});
  end
  combined = vertcat(combined{:});
  if isempty(combined)
    continue
  end
  
  % calculate MIR
  if strcmp(data_type, 'lfp')
    mir_func = @getMirLfp;
  else
    mir_func = @(d) getMirSpike(d, binwidth);
  end
  
  if nargout <= 1
    mir{i} = mir_func(combined);
  else
    mirs = zeros(1,10);
    for j=1:10
      mirs(j) = mir_func(combined);
    end
    mir{i} = mean(mirs);
    se{i} = std(mirs) / sqrt(10);
    if nargout > 2
      quantiles = zeros(1,100);
      for j=1:100
        idxs = randperm(size(combined,1));
        combined.label = combined.label(idxs);
        quantiles(j) = mir_func(combined);
      end
      percentile{i} = quantile(quantiles, args.quantile);
    end
  end
end
if ~strcmp(data_type, 'lfp') || ~args.expand_electrodes
  mir = [mir{:}];
  se = [se{:}];
  percentile = [percentile{:}];
end
if strcmp(data_type, 'spike') && args.average
  mir = mean(mir);
  se = mean(se);
  percentile = mean(percentile);
end
end
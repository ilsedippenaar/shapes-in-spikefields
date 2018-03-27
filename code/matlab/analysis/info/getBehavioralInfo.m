function mir = getBehavioralInfo(dhs, binwidth, delay)
name = dhs(1).monkey_name;
assert(all(strcmp(name, {dhs.monkey_name})));
if strcmp(name, 'zorin')
  h_time = 83; % shape on for 83 ms for zorin
else
  h_time = 120; % shape on for 120 ms for jaws
end

trials = [dhs.trials];
all_vecs = cell(2,numel(trials));
for i=1:numel(trials)
  % eye timepoints
  f = trials(i).sections{3} - delay; % fixate
  s = trials(i).saccade - delay; % saccade
  
  % stimlus timepoints
  h = trials(i).sections{5}; % shape
  n1 = trials(i).sections{4}; % noise onset
  n2 = trials(i).sections{6}; % noise off = trial end?
  % TODO: is n2 really noise off?
  
  % validate trial - make sure times occur in expected places
  if isempty(f) || isempty(s) || ...
      (~isempty(h) && h + h_time > n2)
    continue
  end
  if h
    % align to shape and ensure that h isn't repeated twice in edges
    edges = [fliplr(h:-binwidth:n1-binwidth), h+binwidth:binwidth:n2+binwidth];
    bin_nums = discretize([h,h+h_time], edges);
    stim_vec = zeros(1,numel(edges)-1,'int8');
    stim_vec(bin_nums(1):bin_nums(2)) = 1;
  else
    edges = f:binwidth:n2+binwidth; % align to fixation
    edges = edges(edges > n1-binwidth);
    stim_vec = zeros(1,numel(edges)-1,'int8');
  end
  if s < edges(1) || s >= edges(end) || f > edges(1)
    continue
  end
  eye_vec = zeros(1,numel(edges)-1,'int8');
  eye_vec(discretize(s, edges):end) = 1;
  all_vecs{1,i} = stim_vec;
  all_vecs{2,i} = eye_vec;
end
mir = mutInfo(double([all_vecs{1,:}]), double([all_vecs{2,:}])) / binwidth * 1000;
end
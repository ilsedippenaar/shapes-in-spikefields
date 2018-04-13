function analysisByRow(data, ~, intermKVStore, params)
mt_params = [];
mt_params.tapers = [params.T*params.W, 2*params.T*params.W-1];
mt_params.fpass = [0,params.freq_cutoff];
mt_params.Fs = params.Fs;

num_trials = size(data,1);
out = cell(num_trials,1);
for i=1:num_trials
  out{i} = cohByRow(data(i,:).lfps{1}, data(i,:).electrode_mapping{1}, mt_params); 
end
add(intermKVStore, 'coherence', cat(4,out{:}));
end
%% Parameters
params = [];
% Selecting
params.names = {'shape', 'saccade'};
params.length_postnoise_response = 1001+150;
params.data_selections = [-150,-150;150,150];
% Multitaper
params.T = 0.5;
params.W = 6;
params.window_size = 0.064;
% Plotting
params.freq_cutoff = 120;
%% Initialize
conditions = cell(1,2);
conditions{1}(1).vec = 'shape';
conditions{1}(1).negate = false;
conditions{1}(1).range = [0,1];
conditions{1}(2).vec = 'noise';
conditions{1}(2).negate = true;
conditions{1}(2).range = [-params.length_postnoise_response, 0];

conditions{2} = conditions{1};
conditions{2}(1).vec = 'saccade';

electrode_mapping = cell(96,2);
electrode_mapping(:,1) = num2cell(1:96);

data = cell(1,size(params.data_selections,2));
day_lfps = cell(size(data));
for i=1:numel(data)
  day_lfps{i} = cell(1,numel(dhs));
  for j=1:numel(dhs)
    fprintf('%d / %d\n', j, numel(dhs));
    day_lfps{i}{j} = expandCellArray(dhs(j).electrode_mapping, ...
      dhs(j).getDataSlices('lfp', params.data_selections(:,i), conditions{i}));
  end
  data{i} = combineCellArrays('single', day_lfps{i}{:});
end

combined_indices = cell(size(data));
for i=1:numel(data)
  combined_indices{i} = cell(1,96);
  num_in_day = cellfun(@(c) max(cellfun(@(d) size(d,2), c)), day_lfps{i});
  for j=1:96
    combined_indices{i}{j} = zeros(1,sum(num_in_day));
    curr_idx = 1;
    for k=1:numel(day_lfps{i})
      if ~isempty(day_lfps{i}{k}{j})
        combined_indices{i}{j}(curr_idx:curr_idx+num_in_day(k)-1) = 1;
      end
      curr_idx = curr_idx + num_in_day(k);
    end
    combined_indices{i}{j} = find(combined_indices{i}{j});
  end
end

%% Sample trial trace
%% Traces
for i=1:numel(data)
  lfps = cellArray2mat(data{i},'single');
  plt = plotMeanAndStds(lfps, 'x', params.data_selections(1,i):params.data_selections(2,i)-1);
  plt.Color = 'white';
  xlabel('Time (ms)', 'Interpreter', 'latex');
  ylabel('Voltage ($\mu$V)', 'Interpreter', 'latex');
  saveFigures(plt, fullfile(plot_save_dir, 'poster', sprintf('trace_%s.png', params.names{i})));
end
%% Spectrograms
mt_params = [];
mt_params.tapers = [params.T*params.W, 2*params.T*params.W-1];
mt_params.fpass = [0,params.freq_cutoff];
mt_params.Fs = 1000;
mt_params.trialave = true;
for i=1:numel(data)
  lfps = cellArray2mat(data{i}, 'single');
  [S,t,f] = mtspecgramc(lfps, [params.window_size, params.window_size/2], mt_params);
  plt = plotSpectrogram(S,linspace(params.data_selections(1,i),params.data_selections(2,i),numel(t))/1000,f);
  plt.Color = 'white';
  colormap(plt, 'jet');
  saveFigures(plt, fullfile(plot_save_dir, 'poster', sprintf('spectrogram_%s.png', params.names{i})));
end
%% Coheregrams
mt_params = [];
mt_params.tapers = [params.T*params.W, 2*params.T*params.W-1];
mt_params.fpass = [0,params.freq_cutoff];
mt_params.Fs = 1000;
mt_params.trialave = true;
all_cohgrams = cell(size(data));
for i=1:numel(data)
  C = cohgramc(data{i}{1}(:,1),data{i}{1}(:,1), [params.window_size, params.window_size/2], mt_params);
  cohgram = zeros(size(C));
  all_cohgrams{i} = zeros([size(C), 96, 96]);
  for j=1:96
    fprintf('%d / %d\n', j, 96);
    all_cohgrams{i}(:,:,j,j) = 1;
    for k=j+1:96
      [~,idxs1,idxs2] = intersect(combined_indices{i}{j},combined_indices{i}{k});
      [C,~,~,~,~,t,f] = cohgramc(data{i}{j}(:,idxs1), data{i}{k}(:,idxs2),  [params.window_size, params.window_size/2], mt_params);
      cohgram = cohgram + C;
      all_cohgrams{i}(:,:,k,j) = C;
      all_cohgrams{i}(:,:,j,k) = C;
    end
  end
  cohgram = cohgram / nchoosek(96,2);
  plt = plotSpectrogram(cohgram,linspace(params.data_selections(1,i),params.data_selections(2,i),numel(t))/1000,f,false);
  plt.Color = 'white';
  colormap(plt, 'jet');
  saveFigures(plt, fullfile(plot_save_dir, 'poster', sprintf('coheregram_%s.png', params.names{i})));
end
%% Interelectrode distance coherence
breaks = linspace(0, hypot(9,9)*400, 6);
dist_bins = [breaks(1:end-1);breaks(2:end)];

mt_params = [];
mt_params.tapers = [params.T*params.W, 2*params.T*params.W-1];
mt_params.fpass = [0,params.freq_cutoff];
mt_params.Fs = 1000;
mt_params.trialave = true;
all_cohs = cell(size(data));
for i=1:numel(data)
  C = coherencyc(data{i}{1}(:,1),data{i}{1}(:,1), mt_params);
  cohs = zeros(numel(C),96,96);
  for j=1:96
    fprintf('%d / %d\n', j, 96);
    cohs(:,j,j) = 1;
    for k=j+1:96
      [~,idxs1,idxs2] = intersect(combined_indices{i}{j},combined_indices{i}{k});
      [C,~,~,~,~,f] = coherencyc(data{i}{j}(:,idxs1), data{i}{k}(:,idxs2), mt_params);
      cohs(:,j,k) = C;
      cohs(:,k,j) = C;
    end
  end
  [inter_cohs, num_in_bins, binned_cohs] = calculateInterelecDistCoherence(cohs, dist_bins, electrode_mapping);
  plt = plotInterelecDistCoherence(inter_cohs, f, num_in_bins, dist_bins, 'freq_bin', [0, params.freq_cutoff]);
  plt.Color = 'white';
  saveFigures(plt, fullfile(plot_save_dir, 'poster', sprintf('interelectrode_dist_coh_%s.png', params.names{i})));
  all_cohs{i} = cohs;
end
%% Export coherency for plotting in R
r_data = [];
r_data.binned_cohs = binned_cohs;
r_data.freqs = f;
r_data.dist_bins = dist_bins;
save(fullfile(data_dir, '../R/'), '-struct', r_data);
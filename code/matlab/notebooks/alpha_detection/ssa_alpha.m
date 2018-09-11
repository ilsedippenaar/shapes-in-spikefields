monkey_name = 'jaws';
folder_name = 'alpha_detection/ssa_alpha';

trial_dir = fullfile(data_dir, 'trials', monkey_name);
lfp_dir = fullfile(data_dir, 'lfps', monkey_name);
dh_dir = fullfile(cache_dir, 'data_handlers', monkey_name);
save_dir = fullfile(cache_dir, 'trial_tables', 'all', monkey_name);

if strcmp(monkey_name, 'jaws')
  dates = datetime('2013-05-1'):datetime('2013-07-19');
  learn_start = datetime('2013-06-3');
else
  dates = datetime('2011-7-22'):datetime('2011-12-20');
  learn_start = datetime('2011-11-7');
end
%% 
tbl = makeTrialTable(dates, trial_dir, lfp_dir, dh_dir, save_dir, [], ...
  'clean', true, ...
  'min_reaction_time', params.length_min_reaction_time, ...
  'monkey_name', monkey_name);
tbl = tbl(~strcmp(tbl.result, 'failed'),:);
all_tbl = tbl;
[tbl,t] = select(all_tbl, 'shape');
shape_idx = find(t==0);
%%
tbl.success = strcmp(tbl.result, 'true_positive');
tbl.trace = cellfun(@(c) mean(c,2), tbl.lfps, 'UniformOutput', false);
tbl.lab = grp2idx(tbl.result);
tbl.rxn_time = tbl.saccade - tbl.shape;

ssa = cell(1,size(tbl,1));
contribs = cell(size(ssa));
freqs = cell(size(ssa));
parfor ii=1:size(tbl,1)
  [ssa{ii}, contribs{ii}, freqs{ii}] = getSSAComponents(tbl(ii,:).trace{1}, 1000);
end
tbl.ssa = ssa';
tbl.contribs = contribs';
tbl.ssa_freqs = freqs';
%% 
a_tbl = tbl(contains(tbl.result, {'true_positive', 'false_negative'}),:);
shape_angles = cell(size(a_tbl,1),1);
shape_mag = cell(size(shape_angles));
parfor ii=1:size(a_tbl,1)
  % [~, freq_idx] = min(abs(a_tbl(ii,:).ssa_freqs{1} - 12));
  %shape_idx = a_tbl(ii,:).shape - a_tbl(ii,:).start + 1;
  z = hilbert(a_tbl(ii,:).ssa{1});
  shape_angles{ii} = angle(z(shape_idx,:));
  shape_mag{ii} = a_tbl(ii,:).contribs{1}'.*abs(z(shape_idx,:));
end
%% 
n_bins = 50;
freqs = 1:100;
freq_infos = zeros(1, numel(freqs));
mag_infos = zeros(size(freq_infos));
for i=1:numel(freqs)
  f_idxs = cellfun(@(c) argmin(abs(c - freqs(i))), a_tbl.ssa_freqs);
  freq_angles = zeros(1,numel(f_idxs));
  freq_mags = zeros(size(freq_angles));
  for j=1:numel(f_idxs)
    freq_angles(j) = shape_angles{j}(f_idxs(j));
    freq_mags(j) = shape_mag{j}(f_idxs(j));
  end
  s = histcounts(freq_angles(a_tbl.success), linspace(-pi, pi, n_bins));
  f = histcounts(freq_angles(~a_tbl.success), linspace(-pi, pi, n_bins));
  s = s / sum(s);
  f = f / sum(f);
  [~, freq_infos(i)] = ttest(s,f);
  
  % multiply by -1 if angle would make magnitude 'negative'
  y = freq_mags; % .* (between(wrapTo2Pi(freq_angles), pi/2, 3*pi/2)*-2+1);
  m = min(y);
  M = max(y);
  s = histcounts(y(a_tbl.success), linspace(m, M, n_bins));
  f = histcounts(y(~a_tbl.success), linspace(m, M, n_bins));
  s = s / sum(s);
  f = f / sum(f);
  [~, mag_infos(i)] = ttest(s,f);
end

figure
plot(freqs, freq_infos)
title('Normalized MI of phase angles for S/F vs frequency')
%saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('nmi_phase_%s.png', monkey_name)));

figure
plot(freqs, mag_infos);
title("p val of paired t-test for magnitude vs frequency")
%saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('pval_mag_%s.png', monkey_name)));
% no phase or magnitude relationship at any frequency of SSA components on shape
% detection at time of shape onset
%% Phase effects and reaction time??
n_bins = 10;
freqs = 1:100;
corrs = zeros(1, numel(freqs));
mag_corrs = zeros(size(corrs));
suc_shape_angles = shape_angles(a_tbl.success);
suc_shape_mags = shape_mag(a_tbl.success);
for i=1:numel(freqs)
  f_idxs = cellfun(@(c) argmin(abs(c - freqs(i))), a_tbl(a_tbl.success,:).ssa_freqs);
  freq_angles = zeros(1,numel(f_idxs));
  freq_mags = zeros(size(freq_angles));
  for j=1:numel(f_idxs)
    freq_angles(j) = suc_shape_angles{j}(f_idxs(j));
    freq_mags(j) = suc_shape_mags{j}(f_idxs(j));
  end
  R = corrcoef(freq_angles, a_tbl(a_tbl.success,:).rxn_time);
  corrs(i) = R(2,1);
  
  R = corrcoef(freq_mags, a_tbl(a_tbl.success,:).rxn_time);
  mag_corrs(i) = R(2,1);
end

figure
plot(freqs, corrs)
title("Corr coef of phase and rxn time vs frequency")
%saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('corr_phase_rxn_%s.png', monkey_name)));

figure
plot(freqs, mag_corrs)
title("Corr coef of magnitude and rxn time vs frequency")
%saveFigures(gcf, fullfile(plot_save_dir, folder_name, sprintf('corr_mag_rxn_%s.png', monkey_name)));
% no (linear) correlations of phase angle with reaction time

% approx. 0.15 positive correlation of magnitude of envelope of alpha range
% with reaction time, so higher magnitude ==> slower reaction
% however, this is possibly just due to electrode losing range and Jaws
% getting back over time -- both are caused by time
%%
% figure
% histogram(shape_angles(a_tbl.success));
% hold on
% histogram(shape_angles(~a_tbl.success));
% %%
% dhs = DataHandler.fromDates(dates(1:2), trial_dir, lfp_dir, dh_dir, ...
%   'clean', true, ...
%   'min_reaction_time', params.length_min_reaction_time, ...
%   'monkey_name', monkey_name);
% lfp = mean(dhs(1).lfps,2);
% y = single(lfp(10000:20000));
% y = (y - mean(y)) / std(y);
% [Q,D] = ssacom(y, 100);
% % 3rd reconstructed component has frequency concentration in 5-15 Hz
% z = hilbert(Q(:,1));
% figure,plot(abs(z))
% 
% figure, plot(y)
% hold on
% plot(Q(:,1:3))
% %%

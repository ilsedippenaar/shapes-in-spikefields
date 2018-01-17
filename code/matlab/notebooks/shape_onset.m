%% Initialize
names = {'success', 'failure'};
select_range = [-256, 256];
conds = {
  cell2struct([...
    {'shape', false, [0,1]}; ...
    %{'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]}; ...
    %{'saccade', true, [-1000, params.length_min_reaction_time]}; ... % shapes can come on *after* a saccade too
    {'saccade', false, [params.length_min_reaction_time, 1000]}], {'vec','negate','range'}, 2), ...
  cell2struct([...
    {'shape', false, [0,1]}; ...
    %{'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]}; ...
    {'saccade', true, [-1000, 1000]}], {'vec','negate','range'}, 2)};
freq_bins = 0:params.Fs/diff(select_range):25;

lfps = cell(numel(dhs),numel(conds));
for cond_num=1:numel(conds)
  for i=1:numel(dhs)
    lfps{i,cond_num} = dhs(i).getDataSlices('lfp', select_range, conds{cond_num});
  end
end
% %%
% % TODO: look at saccade - shape histogram to determine min_reaction_time
% success_times = [];
% fail_times = [];
% for dh=dhs
%   success_idxs = contains({dh.trials.result}, "true_positive");
%   fail_idxs = contains({dh.trials.result}, "false_negative");
%   success_times = [success_times cellfun(@(c) c{5}-c{4}, {dh.trials(success_idxs).sections})];
%   fail_times = [fail_times cellfun(@(c) c{5}-c{4}, {dh.trials(fail_idxs).sections})];
% end
% 
% % Phase Distribution with Noise Onset at Phase = 0
% plt = figure('Visible', 'off');
% subplot(1,2,1);
% histogram(mod(success_times, 1000 / 12) * 12 / 1000 * 360, 'BinWidth', 36);
% yl = ylim;
% title('Sucess');
% subplot(1,2,2);
% histogram(mod(fail_times, 1000 / 12) * 12 / 1000 * 360, 'BinWidth', 36);
% ylim(yl);
% title('Failure');
% saveFigures(plt, fullfile(plot_save_dir, 'shape_onset', 'phase_distribution_basic.png'));
% %% Same as above, just with HHT
% % for each day and trial (success/failure), get IMFs for mean LFPs from noise onset to shape stimulus
% % consider IMFs with max frequency response in [10,20]
% % do Hilbert transform on those IMFs and calculate the phase angle at shape onset
% % plot histogram of phase angles for success and failure trials
% 
% % TODO: use amplitdue / power to weight phases by relative importance
% % TODO: low pass filter under 50 Hz
% % TODO: repeat above using SVD / SSA instead of EMD
% 
% var_conds = {{'type', 'lfp', 'trial_result', {'true_positive'}, 'trial_section', {'noise_to_shape'}}, ...
%          {'type', 'lfp', 'trial_result', {'false_negative'}, 'trial_section', {'noise_to_shape'}}};
% freq_range = [10,20];
% var_length_phases = getVariableLengthLfpPhases(var_conds, freq_range);
% %% Plotting variable length
% for cond_num=1:numel(var_conds)
%   plt = figure('Visible', 'off');
%   histogram(wrapTo2Pi(angle(var_length_phases{cond_num}))/pi*180, 12);
%   saveFigures(plt, fullfile(plot_save_dir, 'shape_onset', sprintf('phase_distribution_hilbert_%s.png', names{cond_num})));
% end
%%
phases_at_shape = getStaticLengthLfpPhases(lfps, conds, freq_bins, -select_range(1));
%% Calculate KL divergence for each frequency bin to find the frequencies most correlated with differences in failure / success
kl_divs = zeros(1,numel(freq_bins)-1);
for i=1:numel(freq_bins)-1
  kl_divs(i) = KLDiv(histcounts(wrapTo2Pi(angle(phases_at_shape{1}{i})), 12), ...
                     histcounts(wrapTo2Pi(angle(phases_at_shape{2}{i})), 12));
end
%% Two plots, one for success / failure
for cond_num=1:numel(conds)
  plt = figure('Visible', 'off');
  % 7 corresponds to the 11.7 - 13.6 (high alpha) bin
  histogram(wrapTo2Pi(angle(phases_at_shape{cond_num}{6}))/pi*180, 12);
  xlabel('Phase Angle');
  %title('Success rate for high alpha 11.7-13.6 Hz');
  title('Success rate for high alpha 9.8-11.7 Hz');
  saveFigures(plt, fullfile(plot_save_dir, 'shape_onset', sprintf('phase_distribution_11.7_to_13.6_%s.png', names{cond_num})));
end
%% Sucess rate plots (s/(s+f))
nbins = 8;
success_counts = histcounts(wrapTo2Pi(angle(phases_at_shape{1}{12})), nbins);
failure_counts = histcounts(wrapTo2Pi(angle(phases_at_shape{2}{12})), nbins);
plt = figure('Visible', 'off');
histogram('BinEdges', linspace(0,360,nbins+1), 'BinCounts', success_counts ./ (success_counts + failure_counts))
xlabel('Phase Angle');
ylabel('Success rate');
title('Success rate for high alpha 11.7-13.6 Hz');
saveFigures(plt, fullfile(plot_save_dir, 'shape_onset', 'success_rate_plot_high_alpha.png'), false);
savefig(plt, fullfile(plot_save_dir, 'figs/phase', 'success_rate.fig'));
close(plt);
%%
circ_rtest(success_counts ./ (success_counts + failure_counts) * sum(success_counts+failure_counts)/nbins)
%%
function phases_at_shape = getHHTPhases(x, freq_bins, phase_at)
imfs = emd(x);
imf_freqs = cellfun(@(c) argmax(pwelch(c, hann(512)))/257*1000/2, imfs); % get frequency with greatest power in each IMF
% [~,closest_to_17] = min(abs(imf_freqs - 17));
phases_at_shape = cell(1,numel(freq_bins)-1);
for i=1:numel(freq_bins)-1
  freq_imfs = imfs(and(imf_freqs >= freq_bins(i), imf_freqs < freq_bins(i+1)));
  % imfs  = imfs(closest_to_17);
  phases_at_shape{i} = zeros(1,numel(freq_imfs));
  for j=1:numel(freq_imfs)
    h = hilbert(freq_imfs{j});
    if phase_at == 'end'
      phases_at_shape{i}(j) = h(end);
    else
      phases_at_shape{i}(j) = h(phase_at);
    end
  end
end
end

function phases_at_shape = getVariableLengthLfpPhases(conds, freq_range)
phases_at_shape = cell(1,numel(conds));
for cond_num=1:numel(conds)
  phases_at_shape{cond_num} = cell(1,numel(dhs));
  for i=1:numel(dhs)
    fprintf('Day = %d\n', i);
    day_lfps = dhs(i).select(conds{cond_num}{:});
    lfps = cell(numel(day_lfps{1}),numel(day_lfps));
    for j=1:size(lfps,2) % lfps is num_trials x num_electrodes
      lfps(:,j) = day_lfps{j}';
    end
    valid_trials = cellfun(@numel, lfps(:,1)) > 512;
    lfps = lfps(valid_trials, :);
    phases_at_shape{cond_num}{i} = cell(1,numel(lfps));
    for j=1:size(lfps,1) % for each trial
      lfp = mean(cellArray2mat(lfps(j,:)),2);
      phases_at_shape{cond_num}{i}{j} = getHHTPhases(lfp, freq_range, 'end');
    end
    phases_at_shape{cond_num}{i} = [phases_at_shape{cond_num}{i}{:}];
  end
  phases_at_shape{cond_num} = [phases_at_shape{cond_num}{:}];
end
end

function phases_at_shape = getStaticLengthLfpPhases(lfps, conds, freq_bins, phase_at)
phases_at_shape = cell(1,numel(conds));
for cond_num=1:numel(conds)
  phases_at_shape{cond_num} = cell(1,size(lfps,1));
  for i=1:numel(phases_at_shape{cond_fnum})
    fprintf('Day = %d\n', i);
    day_lfps = lfps{i, cond_num};
    day_lfps = reshape([day_lfps{:}], [size(day_lfps{1}), numel(day_lfps)]); % time x trials x electrode
    day_lfps = mean(day_lfps,3); % time x trials
    phases_at_shape{cond_num}{i} = cell(1,size(day_lfps,2));
    for j=1:size(day_lfps,2) % for each trial
      phases_at_shape{cond_num}{i}{j} = getHHTPhases(day_lfps(:,j), freq_bins, phase_at);
    end
    phases_at_shape{cond_num}{i} = combineCellArrays('complex', phases_at_shape{cond_num}{i}{:});
  end
  phases_at_shape{cond_num} = combineCellArrays('complex', phases_at_shape{cond_num}{:});
end
end

function [pval z] = circ_rtest(alpha, w, d)
%
% [pval, z] = circ_rtest(alpha,w)
%   Computes Rayleigh test for non-uniformity of circular data.
%   H0: the population is uniformly distributed around the circle
%   HA: the populatoin is not distributed uniformly around the circle
%   Assumption: the distribution has maximally one mode and the data is 
%   sampled from a von Mises distribution!
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r, in radians (!)]
%
%   Output:
%     pval  p-value of Rayleigh's test
%     z     value of the z-statistic
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if size(alpha,2) > size(alpha,1)
	alpha = alpha';
end

if nargin < 2
	r =  circ_r(alpha);
  n = length(alpha);
else
  if length(alpha)~=length(w)
    error('Input dimensions do not match.')
  end
  if nargin < 3
    d = 0;
  end
  r =  circ_r(alpha,w(:),d);
  n = sum(w);
end

% compute Rayleigh's R (equ. 27.1)
R = n*r;

% compute Rayleigh's z (equ. 27.2)
z = R^2 / n;

% compute p value using approxation in Zar, p. 617
pval = exp(sqrt(1+4*n+4*(n^2-R^2))-(1+2*n));
end
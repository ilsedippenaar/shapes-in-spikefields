%% Initialize
names = {'success', 'failure'};
select_range = [-256, 0];
conds = {
  cell2struct([...
    {'shape', false, [0,1]}; ...
    {'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]}; ...
    {'saccade', true, [-1000, params.length_min_reaction_time]}; ... % shapes can come on *after* a saccade too
    {'saccade', false, [params.length_min_reaction_time, 1000]}], {'vec','negate','range'}, 2), ...
  cell2struct([...
    {'shape', false, [0,1]}; ...
    {'noise', true, [-params.length_postnoise_response+select_range(1), select_range(2)]}; ...
    {'saccade', true, [-1000, 1000]}], {'vec','negate','range'}, 2)};
freq_bins = 0:params.Fs/diff(select_range):params.Fs/2;
[~,alpha_idx] = min(abs(freq_bins-10));
freq_cutoff = 100;
folder_name = 'jan_grant_app';

lfps = cell(numel(dhs),numel(conds));
for cond_num=1:numel(conds)
  for i=1:numel(dhs)
    lfps{i,cond_num} = dhs(i).getDataSlices('lfp', select_range, conds{cond_num});
  end
end

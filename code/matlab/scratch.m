% a = cell(1,numel(dh.trials));
% for i=1:numel(dh.trials)
%   if ~isempty(dh.trials(i).sections{5}) && ~isempty(dh.trials(i).saccade)
%     a{i} = dh.trials(i).saccade - dh.trials(i).sections{5};
%   end
% end
% a = [a{:}];
% mean(a)
% 
% lfp_slice = dh.getLfpSlices('200,shape,1,0,none,500,500', between(1), between(2));
% plotMeanAndStds(lfp_slice{1})
% hold on
% plot([200 200], [-400 400], 'r--')
%%
% Fs = 40000;
% sample = double(dh.lfps(1:Fs*10,1));
% min_sample = min(sample);
% max_sample = max(sample);
% scaled_sample = (sample-min_sample)*(double(intmax('uint16'))/(max_sample-min_sample)) + double(intmin('int16'));
% wave = audioplayer(int16(sample), Fs);
% play(wave);

%%
% lfp_idxs = 1.45e6:2.7e6;
% tmp = fft(dh.lfps(lfp_idxs,1));
% v = abs(tmp);
% n = numel(v);
% idxs = 1.686e5:1.6875e5;
% idxs = [idxs (n+1-idxs)];
% new_spec = tmp;
% new_spec(~ismember(1:n, idxs)) = 0;
% filtered = real(ifft(new_spec));
% min_filt = min(filtered);
% max_filt = max(filtered);
% scaled_filt = (filtered-min_filt)*(double(intmax('uint16'))/(max_filt-min_filt)) + double(intmin('int16'));
% Fs2 = 4000;
% wave = audioplayer(int16(scaled_filt(1:Fs2*10)), Fs2);
% play(wave);
%%
% lfp_sections = zeros(1000, 89*numel(lfp_slices));
% labels = cell(1,89*numel(lfp_slices));
% for i=1:numel(lfp_slices)
%   for j=1:89
%     lfp_sections(:,(i-1)*89 + j:(i-1)*89+j) = mean(lfp_slices{i}{j},2);
%     labels{(i-1)*89+j} = names{i};
%   end
% end
% lfp_sections = lfp_sections';
% save(fullfile(data_dir,'python','lfp_sections.mat'), 'lfp_sections')
% save(fullfile(data_dir,'python','labels.mat'), 'labels')
% %%
% save(fullfile(data_dir, 'python', 'events.mat'), 'labels', 'all_events')

%%
zz = hilbert(single(lfps{2}(1:10000)));
subplot(5,1,1);
w = diff(unwrap(angle(zz)));
plot(w);
subplot(5,1,2);
plot(angle(zz));
subplot(5,1,3);
tmp = spikes{1}(and(spikes{1} >= start_time, spikes{1} < start_time+10000));
hold on
for i=1:numel(tmp), plot([tmp(i) tmp(i)], [-10 10], 'LineWidth', 0.1, 'Color', 'black'); end
subplot(5,1,4);
tmp2 = zeros(1,10000);
tmp2(tmp-start_time) = 1;
plot(smooth(tmp2,100,'lowess'));
subplot(5,1,5);
plot(abs(zz));
%%
dist_mat = zeros(89,89);
for i=1:89
  dist_mat(i,i) = 1;
  for j=i+1:89
    v1=single(lfps{i});
    v2=single(lfps{j});
    dist_mat(i,j) = dist(v1,v2);
    dist_mat(j,i) = dist_mat(i,j);
  end
end
figure;
imagesc(dist_mat);
colormap hot;
colorbar;
%%
for i=1:get(gcf, 'Number')
  close(gcf);
end
%%
e = [];
for j=1:numel(dhs)
    d = cellfun(@(c) c{5}, {dhs(j).trials.sections}, 'UniformOutput', false);
for i=1:numel(d)
    if ~isempty(d{i}) && ~isempty(dhs(j).trials(i).saccade)
        e = [e, d{i}-dhs(j).trials(i).saccade];
    end
end
end
max(e)
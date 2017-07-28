%% Reaction times (see plots/R)
rxns_and_slopes = regressReactionTimes(dh, [start_time, size(dh.lfps,1)-30*1000]);
csvwrite(fullfile(data_dir, '../R/', 'rxns_and_slopes.csv'), ...
  [repelem(1:size(rxns_and_slopes,3), size(rxns_and_slopes,1))', ...
  reshape(permute(rxns_and_slopes,[1,3,2]), [], 2)]);
%% Cosine distance in time domain
dist = @(v1,v2) (v1'*v2) / (norm(v1)*norm(v2));
dists = zeros(numel(lfp_slices{1}), 2);
for i=1:numel(lfp_slices{1})
  mean_noise_lfp = mean(lfp_slices{contains(names,'noise')}{i},2);
  mean_shape_lfp = mean(lfp_slices{contains(names,'shape')}{i},2);
  mean_saccade_lfp = mean(lfp_slices{contains(names,'saccade')}{i},2);
  dists(i,:) = acos(abs([...
    dist(mean_noise_lfp, mean_shape_lfp), ...
    dist(mean_noise_lfp, mean_saccade_lfp)]));
end
plot(dists(:,1), dists(:,2), '.', 'MarkerSize', 10);
xlim([0,pi/2]);
ylim([0,pi/2]);
hold on;
plot(acos(linspace(0,1,100)),acos(linspace(0,1,100)));
%% Dimension reduction (PCA) in time domain
n = numel(lfp_slices{1});
all_traces = zeros(diff(select_range),3*n);
cond_names = {'noise','shape','saccade'};
for i=1:numel(cond_names)
  all_traces(:,((i-1)*n+1:i*n)) = cellArray2mat(...
    cellfun(@(m) mean(m,2), lfp_slices{contains(names,cond_names{i})}, 'UniformOutput', false));
end
[~,S,V] = svds(all_traces,2);
reduced = V*S;
reduced = reshape(reduced', 2, n, 3);
figure;
hold on;
colors = {'red','green','blue'};
for i=1:numel(cond_names)
  plot(reduced(1,:,i),reduced(2,:,i),'.', 'MarkerSize', 10, 'Color', colors{i}, 'DisplayName', cond_names{i});
end
legend show
%% Calculate Euclidean distance in reduced space
dist = @(v1,v2) hypot(v1(1,:)-v2(1,:), v1(2,:)-v2(2,:));
dists = [...
  dist(reduced(:,:,1), reduced(:,:,2))', ...
  dist(reduced(:,:,1), reduced(:,:,3))'];
plot(dists(:,1),dists(:,2), '.', 'MarkerSize', 0.1)
text(dists(:,1),dists(:,2),arrayfun(@num2str, 1:n,'UniformOutput', false));
xlabel('Distance between noise and shape points');
ylabel('Distance between noise and saccade points');
%% Frequency domain clustering

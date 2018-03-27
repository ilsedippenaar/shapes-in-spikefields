function plt = plotInfoSurface(info_mat, data_type, info_type, bws, delays)
plt = figure('Visible', 'off');
imagesc(flipud(info_mat))
colormap jet
colorbar
xticks(1:2:size(info_mat,2));
xticklabels(delays(1:2:end));
xlabel('Delay (ms)');
yticks(1:2:size(info_mat,1));
yticklabels(bws(end:-2:1));
ylabel('Binwdith (ms)');
info_type(1) = upper(info_type(1));
if ~isempty(data_type)
  data_type(1) = upper(data_type(1));
  title(sprintf('%s %s Information', data_type, info_type));
else
  title(sprintf('%s Information', info_type));
end
end
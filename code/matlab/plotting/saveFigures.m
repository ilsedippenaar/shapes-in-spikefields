function saveFigures(fig_list, filename)
n = numel(fig_list);
if n == 1
  [~, ~, type] = fileparts(filename);
  export_fig(fig_list(1), filename, sprintf('-%s', type(2:end)), '-a1');
  close(fig_list(1));
  return
end

filenames = cell(1, n);
for i=1:n
  filenames{i} = [tempname '.pdf'];
  export_fig(fig_list(i), filenames{i}, '-pdf', '-a1');
end
if exist(filename, 'file') == 2
  delete(filename);
end
append_pdfs(filename, filenames{:});
cellfun(@(f) delete(f), filenames);
for i=1:n
  close(gcf);
end
end
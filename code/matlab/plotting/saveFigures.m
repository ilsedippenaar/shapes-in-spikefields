function saveFigures(fig_list, filename, close_figs)
%SAVEFIGURES Saves a figure (or figures) to a file
% The parameter fig_list may contain 1 or many MATLAB figures. The default
% output format is pdf, but can be changed by specifying the extension as
% part of filename (for example, "plot1.png"). However, only pdfs will be
% exported if fig_list contains more than one figure.
% 
% close_figs: boolean, default true
if nargin < 3
  close_figs = true;
end
n = numel(fig_list);
if n == 1
  [~, ~, type] = fileparts(filename);
  export_fig(fig_list(1), filename, sprintf('-%s', type(2:end)), '-a1');
  if close_figs
    close(fig_list(1));
  end
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
if close_figs
  for i=1:n
    close(gcf);
  end
end
end
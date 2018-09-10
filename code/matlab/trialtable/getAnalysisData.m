function [lfps,spikes] = getAnalysisData(tbl, analysis_window, align_to)
  idxs = analysis_window(1):analysis_window(2)-1;
  analysis_lfps = rowfun(@(x, start, align) {x{1}(idxs+align-start+1,:)}, ...
    tbl, 'InputVariables', {'lfps', 'start', align_to});
  lfps = table2cell(analysis_lfps);

  select_spikes = @(x, align, range) x(binarySearch(x, range(1), ']', true):binarySearch(x, range(2), '(', true))-align;
  analysis_spikes = rowfun(@(x, align) funcInCells(x, select_spikes, [], {align,idxs([1,end])+align}), ...
    tbl, 'InputVariables', {'spikes', align_to});
  spikes = table2cell(analysis_spikes);
end
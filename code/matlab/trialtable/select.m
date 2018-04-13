function out = select(tbl, type, varargin)
if isnumeric(type)
  out = tbl(type, :);
elseif strcmp(type, 'shape')
  p = inputParser;
  p.addParameter('analysis_window', [-256,256]);
  p.addParameter('min_reaction_time', 200);
  p.addParameter('postnoise_response', 260);
  p.parse(varargin{:});
  args = p.Results;
  valid_trials = [~isnan(tbl.shape), ...
    between(tbl.shape + args.analysis_window, tbl.start, tbl.stop) ... 
    ~between(tbl.noise, tbl.shape - args.postnoise_response + args.analysis_window(1), ...
                        tbl.shape + args.analysis_window(2)), ...
    ~between(tbl.saccade , tbl.shape - 1000, ...
                           tbl.shape + args.min_reaction_time + 1)]; % +1 to agree with DataHandler's definition
  valid_trials = all(valid_trials,2);
  out = tbl(valid_trials,:);
  [out.lfps, out.spikes] = getAnalysisData(out, args.analysis_window, 'shape');
elseif strcmp(type, 'noise')
  p = inputParser;
  p.addParameter('analysis_window', [0,256]+260);
  p.parse(varargin{:});
  args = p.Results;
  valid_trials = [~isnan(tbl.noise), ...
    between([tbl.noise, tbl.noise + args.analysis_window(2)], ...
              tbl.start, tbl.stop)];
  valid_trials = all(valid_trials,2);
  out = tbl(valid_trials,:);
  [out.lfps, out.spikes] = getAnalysisData(out, args.analysis_window, 'noise');
else
  out = [];
end

  function [lfps,spikes] = getAnalysisData(t, analysis_window, align_to)
    idxs = analysis_window(1):analysis_window(2)-1;
    analysis_lfps = rowfun(@(x, start, align) {x{1}(idxs+align-start+1,:)}, ...
      t, 'InputVariables', {'lfps', 'start', align_to});
    lfps = table2cell(analysis_lfps);

    select_spikes = @(x, align, range) x(binarySearch(x, range(1), ']', true):binarySearch(x, range(2), '(', true))-align;
    analysis_spikes = rowfun(@(x, align) funcInCells(x, select_spikes, [], {align,idxs([1,end])+align}), ...
      t, 'InputVariables', {'spikes', align_to});
    spikes = table2cell(analysis_spikes);
  end
end
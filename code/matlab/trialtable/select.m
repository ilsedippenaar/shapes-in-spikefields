function [out,t] = select(tbl, type, varargin)
if isnumeric(type)
  out = tbl(type, :);
elseif strcmp(type, 'preshape')
  p = inputParser;
  p.addParameter('analysis_window', [-256,0]);
  p.addParameter('postnoise_response', 0);
  p.parse(varargin{:});
  args = p.Results;
  valid_trials = [~isnan(tbl.shape), ...
    between(tbl.shape + args.analysis_window, tbl.start, tbl.stop), ... 
		tbl.noise < tbl.shape - args.postnoise_response + args.analysis_window(1)];
  valid_trials = all(valid_trials,2);
  out = tbl(valid_trials,:);
  [out.lfps, out.spikes] = getAnalysisData(out, args.analysis_window, 'shape');
  t = analysis_window(1):analysis_window(2)-1;
elseif strcmp(type, 'prenoise')
  p = inputParser;
  % 500 ms because the vast majority of noises come on < 512 ms after
  % fixation
  p.addParameter('analysis_window', [-500,0]);
  p.parse(varargin{:});
  args = p.Results;
  valid_trials = [~isnan(tbl.noise), ...
    between(tbl.noise + args.analysis_window, tbl.start, tbl.stop), ...
    tbl.noise + args.analysis_window(1) > tbl.fixate];
  valid_trials = all(valid_trials,2);
  out = tbl(valid_trials,:);
  [out.lfps, out.spikes] = getAnalysisData(out, args.analysis_window, 'noise');
  t = analysis_window(1):analysis_window(2)-1;
elseif strcmp(type, 'postnoise')
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
  t = analysis_window(1):analysis_window(2)-1;
elseif strcmp(type, 'shape')
  p = inputParser;
  p.addParameter('analysis_window', [-400,100]);
  p.addParameter('postnoise_response', 0);
  p.parse(varargin{:});
  args = p.Results;
  valid_trials = [~isnan(tbl.shape), ...
    between(tbl.shape + args.analysis_window, tbl.start, tbl.stop), ... 
		tbl.noise < tbl.shape - args.postnoise_response + args.analysis_window(1), ...
    ~between([tbl.saccade, tbl.noise], tbl.shape + args.analysis_window(1), ...
                          tbl.shape + args.analysis_window(2))];
  valid_trials = all(valid_trials,2);
  out = tbl(valid_trials,:);
  [out.lfps, out.spikes] = getAnalysisData(out, args.analysis_window, 'shape');
  t = analysis_window(1):analysis_window(2)-1;
elseif strcmp(type, 'saccade')
  p = inputParser;
  p.addParameter('analysis_window', [-300,50]);
  p.parse(varargin{:});
  args = p.Results;
  valid_trials = [~isnan(tbl.saccade), ...
    tbl.noise < tbl.saccade + args.analysis_window(1), ...
    between(tbl.saccade + args.analysis_window, ...
              tbl.start, tbl.stop)];
  valid_trials = all(valid_trials,2);
  out = tbl(valid_trials,:);
  [out.lfps, out.spikes] = getAnalysisData(out, args.analysis_window, 'saccade');
  t = analysis_window(1):analysis_window(2)-1;
else
  out = [];
end
end

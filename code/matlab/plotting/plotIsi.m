function plt=plotIsi(spikes, varargin)
p = inputParser;
p.addParameter('log_transform', true, @islogical);
p.addParameter('nbins', 100);
p.addParameter('size', [10 10], @(x) numel(x) == 2);
p.addParameter('spaces_to_exclude', [1 10 91 100]);
p.addParameter('subplot_height', 5000);
p.parse(varargin{:});
args = p.Results;

curr_idx = 1;
plt = figure('Visible', 'off', ...
  'Units', 'points', ...
  'Position', [0 0 args.size(1)*150 args.size(2)*150]);

n = args.size(1)*args.size(2);
for i=1:n
  subplot(args.size(1), args.size(2), i);
  ylim([0 args.subplot_height]);
  if any(args.spaces_to_exclude == i)
    rectangle('FaceColor', 'black');
  else
    isi = diff(spikes{curr_idx});
    if args.log_transform
      histogram(log(isi), args.nbins);
    else
      histogram(isi, args.nbins);
    end
    curr_idx = curr_idx + 1;
  end
  set(gca, 'Visible', 'off');
end
end
function [is_diff, p] = circDistribTest(varargin)
% from N. I. Fisher, Statistical analysis of circular data, 1995, pg 122
ni = cellfun(@numel, varargin)';
assert(all(ni > 10)); % make sure sample size is large enough

for i=1:numel(varargin)
  if isreal(varargin{i})
    varargin{i} = exp(1i*varargin{i});
  end
  varargin{i} = varargin{i}(:);
end
all_vals = vertcat(varargin{:});
circ_ranks = tiedrank(angle(all_vals))*2*pi/sum(ni);
idxs = repUneven(1:numel(varargin), ni)';
all_vals = abs(all_vals) .* exp(1i*circ_ranks);
M = accumarray(idxs, all_vals, [], @sum);
W = 2*sum((abs(M) .^ 2) ./ ni);

% [C,S] = deal(zeros(1,numel(varargin)));
% varargin = cellfun(@(c) mod(c(:), 2*pi), varargin, 'UniformOutput', false);
% all_vals = vertcat(varargin{:});
% ranks = tiedrank(all_vals);
% idxs = repUneven(1:numel(varargin), ni)';
% circ_ranks = accumarray(idxs, ranks*2*pi/sum(ni), [], @(x) {x}); 
% for i=1:numel(circ_ranks)
%   C(i) = sum(cos(circ_ranks{i}));
%   S(i) = sum(sin(circ_ranks{i}));
% end
% W = 2*sum((C.^2 + S.^2) ./ ni);
crit_val = chi2inv(0.95, numel(varargin));
is_diff = W > crit_val;
p = chi2cdf(W, numel(varargin), 'upper');
end
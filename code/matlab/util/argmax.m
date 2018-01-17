function idxs = argmax(x)
% This is again due to matlab silliness. This allows argmax functionality
% inside of anonymous functions
[~,idxs] = max(x);
end
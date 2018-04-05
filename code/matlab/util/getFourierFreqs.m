function f = getFourierFreqs(v, Fs, axis)
if nargin < 3
  axis = find(size(v) ~= 1, 1); % first non-singleton dimension
end
n = size(v, axis);
f = 0:Fs/n:Fs/2;
end
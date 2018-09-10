function [c, contribs, freqs] = getSSAComponentsDecomp(v, Fs)
% algorithm from Bonizzi et al, IEEE EMBS 2012
contrib_left = 1; % the rest of the original signal that's left
y = v;
c = {};
contribs = 1;
while contribs(end) > 0.0001 % 0.01 percent
  y = (y - mean(y)) / std(y); % zero mean and normalize
  [~, idx] = max(abs(fft(y)));
  if idx <= 1
    L = floor(Fs / 2);
  else
    L = floor(numel(y) / idx);
  end
  [Q,D] = ssacom(y, L);
  D = D / sum(D);
  
  c{end+1} = Q(:,1); %#ok
  contribs(end+1) = D(1)*contrib_left; %#ok
  contrib_left = contrib_left * (1-D(1));
  y = y - Q(:,1);
end
contribs = contribs(2:end);
freqs = (cellfun(@(c) argmax(abs(fft(c))), c)-1)*Fs/numel(y);
c = [c{:}];
end
function y = getFreqSuccessBootstrap(s_f, labs)
s = countcats(categorical(labs(s_f,:)));
f = countcats(categorical(labs(~s_f,:)));
y = s ./ (s + f);
end
a=emd(single(dhs(1).lfps(1:100000,1)), 13);
a=cellfun(@(c) c', a, 'UniformOutput', false);

figure;
for i=1:15
    subplot(8,2,i);
    [p,f] = pwelch(a{i}, hann(4096), 2048, 4096, 1000);
    plot(f(f<100), p(f < 100));
    set(gca, 'YTickLabel', []);
end
close all; clear

x = [1:6 5 4:20 33 19:3 4:12 11:0 0:11 11 11 11:19 18:1]';
wlen = 5;
nstds = 3;

[peak_idxs, xmean, xstd] = peak_detect(x, wlen, nstds);

plot(x); hold all
plot(peak_idxs, x(peak_idxs),'o')

x = [1:6 5 4:20 19:3 4:12 11:0 0:11 11 11 11:19 18:1]';
wlen = 5;
nstds = 3;
fs = 20e6;

[tdoas, corr_mag_sq, peak_idxs, lags] = estimate_delay(x, wlen, nstds, fs);

plot(corr_mag_sq); hold all
plot(peak_idxs, corr_mag_sq(peak_idxs),'o')
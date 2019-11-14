function [tdoas, corr_mag_sq, peak_idxs, lags] = ...
    estimate_delay(x, wlen, nstds, fs, show_plots)

numrefs = size(x,2);
Ts = 1/fs;

% Cross correlate all receivers against the first
for ii = 2:numrefs
    [corr_out(:,ii-1), lags(:,ii-1)] = xcorr(x(:,ii), x(:,1));
end

% We consider the magnitudue squared of xcorr for peak detection
corr_mag_sq = abs(corr_out).^2;

% Find the peaks within the window of length wlen
peak_idxs = peak_detect(corr_mag_sq, wlen, nstds, show_plots);

tdoas = lags(peak_idxs)*Ts;
lags = lags(peak_idxs);

% plot(corr_out)
% legend('12','13')


end


function [tdoas, corr_mag_sq, peak_idxs, lags, lags_full, num_samps_from_peak] = ...
    get_tdoa(x, wlen, nstds, fs, percent_of_peak, show_plots)

[nsamps, numrefs] = size(x);
Ts = 1/fs;

% Cross correlate all receivers against the first
corr_out = zeros(2*nsamps-1, numrefs-1);
xmeans = mean(x);
x = x - xmeans;
for ii = 2:numrefs
    [corr_out(:,ii-1), lags(:,ii-1)] = xcorr(x(:,ii), x(:,1));
end

% We consider the magnitudue squared of xcorr for peak detection
corr_mag_sq = abs(corr_out).^2;

% Find the peaks within the window of length wlen
[peak_idxs, num_samps_from_peak] = ...
    peak_detect(corr_mag_sq, wlen, nstds, percent_of_peak, show_plots);

if sum(isnan(peak_idxs)) > 0 
    nanidx = isnan(peak_idxs);
    goodidx = ~isnan(peak_idxs);
%     fprintf(1,'\n\nNot enough peaks detected\n\n')
    tdoas(nanidx) = NaN;
    tdoas(goodidx) = lags(goodidx)*Ts;
    lags(nanidx) = NaN;
    lags(goodidx) = lags(goodidx);
else
    tdoas = lags(peak_idxs)*Ts;
    lags_full = lags; % debug sinc interp
    lags = lags(peak_idxs);
end


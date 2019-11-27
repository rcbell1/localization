function [tdoas, corr_mag_sq, peak_idxs, lags] = ...
    get_tdoa(x, wlen, nstds, fs, show_plots)

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

if sum(isnan(peak_idxs)) > 0 
    nanidx = isnan(peak_idxs);
    goodidx = ~isnan(peak_idxs);
    fprintf(1,'\n\nNot enough peaks detected\n\n')
    tdoas(nanidx) = NaN;
    tdoas(goodidx) = lags(goodidx)*Ts;
    lags(nanidx) = NaN;
    lags(goodidx) = lags(goodidx);
else
    tdoas = lags(peak_idxs)*Ts;
    lags = lags(peak_idxs);
end



% plot(corr_out)
% legend('12','13')


end


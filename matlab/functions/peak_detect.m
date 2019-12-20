function [peak_idxs, num_samps_from_peak, xmean, xstd] = ...
    peak_detect(corr_mag_sq, wlen, nstds, percent_of_peak, show_plots)
% this is for positive valued series only, as is the case after a magnitude
% squared operation. wlen should be odd

% xm = movmax(corr_mag_sq,wlen);    % smooths the series to reject false peaks
% [~,idx_max] = max(xm);  % index where the max inside window begins
[nsamps, npairs] = size(corr_mag_sq);
[xmax, idx_max] = max(corr_mag_sq);  % index where the max inside window begins
xmax = max(xmax.');
% idx_max = idx_max + floor(wlen/2);  % the center of the window is true max

% [val, idxs] = max(x); % simple but will result in false peaks when a true
                        % peak does not exist

% check that the peak is truely a global outlier peak and not local
xmean = mean(corr_mag_sq);
xstd = std(corr_mag_sq);
norm_corr_mag_sq = corr_mag_sq - xmean;

peak_idxs = [];
for ii = 1:npairs
%     if (corr_mag_sq(idx_max(ii),ii) - xmean(ii)) > nstds*xstd(ii)
    if (norm_corr_mag_sq(idx_max(ii),ii) - xmean(ii)) > nstds*xstd(ii)
        peak_idxs(ii) = idx_max(ii);
    else
        peak_idxs(ii) = NaN;
    end
end

% Find the number of samples away from peaks we need to go for the peak 
% value to drop by percent_of_peak
for nn = 1:npairs
    if isnan(peak_idxs(nn))
        num_samps(nn) = nan;
    else
        kk = 1;
        while corr_mag_sq(peak_idxs(nn)+kk,nn) > ...
                corr_mag_sq(peak_idxs(nn),nn)*percent_of_peak
            kk = kk + 1;
        end
        num_samps(nn) = kk;
    end
end
num_samps_from_peak = max(num_samps);

if show_plots == 1
    figure
    plot(norm_corr_mag_sq,'.-'); hold all
    xmin = inf;
    for ii = 1:npairs
        plot([1 length(corr_mag_sq)], [nstds*xstd(ii) nstds*xstd(ii)], '--');
        if ~isnan(peak_idxs(ii))
            plot(peak_idxs(ii), norm_corr_mag_sq(peak_idxs(ii),ii), 'o')
            if xmin > min(norm_corr_mag_sq(peak_idxs(ii)-1,ii),norm_corr_mag_sq(peak_idxs(ii)+1,ii))
                xmin = min(norm_corr_mag_sq(peak_idxs(ii)-1,ii),norm_corr_mag_sq(peak_idxs(ii)+1,ii));
            end
        else
            xmin = 0;
        end
    end
    title('Cross Correlation Peak Location Identification')
    xlabel('Sample Number')
    ylabel('|Correlator Out|^2')
    if ~isnan(num_samps_from_peak)
        axis([-inf inf 0 1.1*xmax])
    else
        axis([-inf inf 0 2.1*nstds*max(xstd.')])
    end

    
    axes('position',[0.6 0.5 0.27 0.3])
    plot(norm_corr_mag_sq,'.-'); hold all
    idx_left = inf;
    idx_right = -inf;
    for ii = 1:npairs
        if ~isnan(peak_idxs(ii))
            plot(peak_idxs(ii), norm_corr_mag_sq(peak_idxs(ii),ii), 'o')
            if idx_left > peak_idxs(ii)
                idx_left = peak_idxs(ii);
            end
            if idx_right < peak_idxs(ii)
                idx_right = peak_idxs(ii);
            end
        else
            idx_left = 0;
            idx_right = nsamps;
        end
    end
    if ~isnan(num_samps_from_peak)
        axis([idx_left-round(num_samps_from_peak/2) idx_right+round(num_samps_from_peak/2) 0.95*xmin 1.05*xmax])
    else
        axis([idx_left idx_right 0.95*xmin 1.05*xmax])
    end
    title('Zoom Peaks')
    xlabel('Sample Number')
    ylabel('|Correlator Out|^2')
end

end


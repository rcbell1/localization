function [peak_idxs, los_peak_lags, num_samps_from_peak, xmean, xstd] = ...
    peak_detect(corr_mag, lags, percent_of_peak, show_plots)
% this is for positive valued series only, as is the case after a magnitude
% squared operation. wlen should be odd

[nsamps, npairs] = size(corr_mag);

% My original way of peak detection, very basic
% [xmax, idx_max] = max(corr_mag);  % index where the max inside window begins
% xmax = max(xmax.');
% 
% % check that the peak is truely a global outlier peak and not local
% xmean = mean(corr_mag);
% xstd = std(corr_mag);
% norm_corr_mag = corr_mag - xmean;
% 
% peak_idxs = [];
% for ii = 1:npairs
% %     if (corr_mag(idx_max(ii),ii) - xmean(ii)) > nstds*xstd(ii)
%     if (norm_corr_mag(idx_max(ii),ii) - xmean(ii)) > nstds*xstd(ii)
%         peak_idxs(ii) = idx_max(ii);
%     else
%         peak_idxs(ii) = NaN;
%     end
% end

% Using MATLABs builtin peak finder method
corr_mag_norm = corr_mag./max(corr_mag); % normalize peak heights
peak_thresh = 0.2;
peak_chooser = 1; % 0 chooses peak with smallest TDOA, 1 chooses highest peak
for ii = 1:npairs
    if peak_chooser == 0
        [peak_vals, peak_idxs2] = findpeaks(corr_mag_norm(:,ii), ...
            'MinPeakHeight', peak_thresh);
        [~,los_peak_idx(ii)] = min(abs(lags(peak_idxs2,1)));
        los_peak_lags(ii) = lags(peak_idxs2(los_peak_idx(ii)),1);
        selected_peaks(ii) = peak_idxs2(los_peak_idx(ii));
    else
        [peak_vals,los_peak_idx(ii)] = max(corr_mag_norm(:,ii));
        peak_idxs2 = los_peak_idx(ii);
        los_peak_lags(ii) = lags(los_peak_idx(ii));
        selected_peaks(ii) = los_peak_idx(ii);
    end

        peak_vals_sv{ii} = peak_vals;
        peak_idxs2_sv{ii} = peak_idxs2;
end
% findpeaks(corr_mag_norm(:,1), 'MinPeakHeight', 0.3, 'Annotate', 'extents')

% figure
% subplot(2,1,1)
% plot(lags(:,1),corr_mag_norm(:,1)); hold on
% plot(lags(peak_idxs2_sv{1},1), peak_vals_sv{1},'vb')
% % plot(lags(peak_idxs2_sv{1}(los_peak_idx(1)),1),peak_vals_sv{1}((los_peak_idx(1))),'rv','markerfacecolor','b')
% % plot([lags(1) lags(end)], [peak_thresh peak_thresh], 'k--')
% axis([-inf inf 0 1.2])
% title('Pair 12')
% xlabel('Lag (Samples)')
% subplot(2,1,2)
% plot(lags(:,2),corr_mag_norm(:,2)); hold on
% plot(lags(peak_idxs2_sv{2},2), peak_vals_sv{2},'vb')
% % plot(lags(peak_idxs2_sv{2}(los_peak_idx(2)),2),peak_vals_sv{2}((los_peak_idx(2))),'rv','markerfacecolor','b')
% % plot([lags(1) lags(end)], [peak_thresh peak_thresh], 'k--')
% axis([-inf inf 0 1.2])
% title('Pair 13')
% xlabel('Lag (Samples)')

peak_idxs = selected_peaks;

% % this is just for debug
% for ii = 1:npairs
%     if peak_idxs(ii) ~= selected_peaks(ii)
%         stop_here = 1;
%     end
% end


% Find the number of samples away from peaks we need to go for the peak 
% value to drop by percent_of_peak
for nn = 1:npairs
    if isnan(peak_idxs(nn))
        num_samps(nn) = nan;
    else
        kk = 1;
        while corr_mag(peak_idxs(nn)+kk,nn) > ...
                corr_mag(peak_idxs(nn),nn)*percent_of_peak
            kk = kk + 1;
        end
        num_samps(nn) = kk;
    end
end
num_samps_from_peak = max(num_samps);

if show_plots == 1
    figure
    plot(norm_corr_mag,'.-'); hold all
    xmin = inf;
    for ii = 1:npairs
        plot([1 length(corr_mag)], [xstd(ii) xstd(ii)], '--');
        if ~isnan(peak_idxs(ii))
            plot(peak_idxs(ii), norm_corr_mag(peak_idxs(ii),ii), 'o')
            if xmin > min(norm_corr_mag(peak_idxs(ii)-1,ii),norm_corr_mag(peak_idxs(ii)+1,ii))
                xmin = min(norm_corr_mag(peak_idxs(ii)-1,ii),norm_corr_mag(peak_idxs(ii)+1,ii));
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
        axis([-inf inf 0 2.1*max(xstd.')])
    end

    
    axes('position',[0.6 0.5 0.27 0.3])
    plot(norm_corr_mag,'.-'); hold all
    idx_left = inf;
    idx_right = -inf;
    for ii = 1:npairs
        if ~isnan(peak_idxs(ii))
            plot(peak_idxs(ii), norm_corr_mag(peak_idxs(ii),ii), 'o')
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


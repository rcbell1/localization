function [coords, bias_coords, covar_coords, mse_coords, tdoas_coarse, ...
    tdoas_refined, prob_correlation, prob_detection, Ntrials, unique] = ...
    get_multiple_single_emitters_fromfile(file_path, targetPos, refPos, ...
    wlen, nstds, percent_of_peak, apply_calibration, calib_path, show_plots)

[numdims, numrefs] = size(refPos);
numpairs = numrefs - 1;         % unique pairs of receivers

% Get the signal emitted by the target
[yblock, fs, bounds] = load_signals_fromfile(file_path, show_plots);
Npulses = length(bounds)-1;
% Npulses = 20;
Ntrials = Npulses;

avg_coords = [0;0];
MSE_coords = [0 0;0 0];
unique_avg = 0;
total_missed_peaks = 0;
detection_count = 0;
Ntrials_with_peak = Ntrials;    % if a peak is missed we skip a trial
coords = nan(numdims, Ntrials);
tdoas_coarse = nan(Ntrials, numpairs);
tdoas_refined = nan(Ntrials, numpairs);
for nn = 1:Ntrials
    
    y = yblock(bounds(nn):bounds(nn+1)-1,:);
        
    % Estimate the delay using the received signals
    [tdoas_coarse(nn,:), corr_mag_sq, peak_idxs, lags, lags_full, num_samps_from_peak] = ...
        get_tdoa(y, wlen, nstds, fs, percent_of_peak, apply_calibration, calib_path, show_plots);

    if sum(isnan(peak_idxs)) ~= numpairs
        detection_count = detection_count + 1;
    end
    peak_count = sum(~isnan(peak_idxs));
    missed_peak_count = sum(isnan(peak_idxs)); 
    total_missed_peaks = total_missed_peaks + missed_peak_count;
    
    if missed_peak_count == 0
        % Refine the tdoas using a super resolution algorithm
        numpairs = size(corr_mag_sq,2); % number of correlated pairs of receivers  
        lags_new = zeros(2*num_samps_from_peak+1, numpairs);
        corr_peak_samples = zeros(2*num_samps_from_peak+1,numpairs);
        for ii = 1:numpairs
            sample_idxs = (peak_idxs(ii)-num_samps_from_peak:peak_idxs(ii)+num_samps_from_peak).';
            lags_new(:,ii) = (lags(ii)-num_samps_from_peak:lags(ii)+num_samps_from_peak).';
            if isnan(sample_idxs)
                corr_peak_samples(:,ii) = NaN;
            else
                corr_peak_samples(:,ii) = corr_mag_sq(sample_idxs, ii);
            end
        end
        corr_peak_samples = {corr_peak_samples; corr_mag_sq};    % temporary while debugging sinc interp
        lags_new = {lags_new; lags_full};
        [tdoas_refined(nn,:), diffs] = refine_tdoa(lags_new, corr_peak_samples, fs, show_plots);

        % Feed the refined TDOAs to a localization algorithm
        [coords(:,nn), unique] = geo_lsq(refPos, tdoas_refined(nn,:));
%         [coords(:,nn), unique] = geo_sphere_int(refPos, tdoas_refined);

        % Compute statistical performance metrics
        avg_coords = avg_coords + coords(:,nn); 
        MSE_coords = MSE_coords + (targetPos - coords(:,nn))*(targetPos-coords(:,nn))';
        unique_avg = unique_avg + unique;
        
    else
        Ntrials_with_peak = Ntrials_with_peak - 1;
        tdoas_refined(nn,:) = nan*ones(1,numpairs);
    end
end
avg_coords = avg_coords/Ntrials_with_peak;
bias_coords = avg_coords - targetPos;
MSE_coords = MSE_coords/Ntrials_with_peak;
covar_coords = MSE_coords - bias_coords*bias_coords';
mse_coords = trace(MSE_coords);
unique_avg = unique_avg/Ntrials_with_peak;
unique = unique_avg;

prob_correlation = 1 - total_missed_peaks/(numpairs*Ntrials);
prob_detection = detection_count/Ntrials;
end


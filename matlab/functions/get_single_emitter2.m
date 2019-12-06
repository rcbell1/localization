function [coords, bias_coords, covar_coords, mse_coords, tdoas, ...
    prob_correlation, prob_detection, unique] = get_single_emitter2(...
    targetPos, refPos, Ntrials, tx_pwr_dbm, fc, fs, fsym, Nsym, span, ...
    sps, beta, wlen, nstds, percent_of_peak, show_plots)

[numdims, numrefs] = size(refPos);
numpairs = numrefs - 1;         % unique pairs of receivers

% Generate the signal emitted by the target
[x, noise_bw] = generate_signal2(Nsym, fsym, sps, span, beta, show_plots); 

%%%%%%% NEED TO RETHINK THE SAMPLE RATE CHANGE HERE. IF SAMPLE RATE IS SET
%%%%%%% LOWER THAN SYMBOL RATE A RESAMPLE NEEDS TO HAPPEN. CURRENTLY THE
%%%%%%% SIGNAL SAMPLE RATE IS SPS*FSYM. THIS SHOULD BE RESAMPLED TO MATCH
%%%%%%% WHATEVER THE RX SAMPLE RATE IS SET TO
[P,Q] = rat(fs/(fsym*sps));
y1 = resample(x, P, Q);

% Add proper delays that correspond to target and emitter locations
[y2, true_tdoas, ranges] = add_delay2(y1, targetPos, refPos, ...
    fs, show_plots);

% y2=y1;

avg_coords = [0;0];
MSE_coords = [0 0;0 0];
unique_avg = 0;
total_missed_peaks = 0;
detection_count = 0;
Ntrials_with_peak = Ntrials;    % if a peak is missed we skip a trial
coords = [inf;inf];
for nn = 1:Ntrials
%         if sum(sum(isnan(MSE_coords)))
%             t = 1;
%         elseif sum(sum(isnan(avg_coords)))
%             t = 2;
%         elseif sum(sum(isnan(coords)))
%             t = 3;
%         end
    % Add noise at the proper SNR levels for free space path losses
    y3 = add_noise(y2, tx_pwr_dbm, noise_bw, fc, ranges, show_plots);

    % Estimate the delay using the received signals
    [tdoas, corr_mag_sq, peak_idxs, lags, num_samps_from_peak] = ...
        get_tdoa(y3, wlen, nstds, fs, percent_of_peak, show_plots);

    %%%% NEED TO TAKE CARE OF MISSING PEAKS GRACEFULLY HERE. SEND BACK A
    %%%% SIGNAL AND PLOT A NO DETECTION COLOR ON THE HEAT MAP. PROBABLY
    %%%% MAKES SENSE TO BYPASS EVERYTHING AFTER THIS POINT WHEN A PEAK
    %%%% DETECT IS MISSING
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
        [tdoas2, diffs] = refine_tdoa(lags_new, corr_peak_samples, fs, show_plots);

        % Feed the refined TDOAs to a localization algorithm
    %     [coords, unique] = geo_lsq(refPos, tdoas2);
        [coords, unique] = geo_sphere_int(refPos, tdoas2);

        % Compute statistical performance metrics
        avg_coords = avg_coords + coords; 
        MSE_coords = MSE_coords + (targetPos - coords)*(targetPos-coords)';
        unique_avg = unique_avg + unique;
        
%         if sum(sum(isnan(MSE_coords)))
%             t = 1;
%         elseif sum(sum(isnan(avg_coords)))
%             t = 2;
%         elseif sum(sum(isnan(coords)))
%             t = 3;
%         end
    else
        Ntrials_with_peak = Ntrials_with_peak - 1;
        tdoas2 = nan*ones(1,numpairs);
%         coords = nan*ones(numdims,1);
    end
end
avg_coords = avg_coords/Ntrials_with_peak;
bias_coords = avg_coords - targetPos;
MSE_coords = MSE_coords/Ntrials_with_peak;
covar_coords = MSE_coords - bias_coords*bias_coords';
mse_coords = trace(MSE_coords);
unique_avg = unique_avg/Ntrials_with_peak;
unique = unique_avg;



tdoas = tdoas2;
prob_correlation = 1 - total_missed_peaks/(numpairs*Ntrials);
prob_detection = detection_count/Ntrials;
end


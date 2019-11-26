function [coords, bias_coords, covar_coords, mse_coords, tdoas, unique] = ...
    get_single_emitter(targetPos, refPos, Ntrials, tx_pwr_dbm, fc, fs, ...
    fsym, Nsym, span, sps, beta, fhigh, wlen, nstds, show_plots)

% fratio = ceil(fhigh/fs);
% sps_high = sps*fratio;      % add delays at high time resolution
% sps_high = fhigh/(fs/sps);  % samples per symbol at high rate
sps_high = fhigh/fsym;      % samples per symbol at high rate
resample_rate = fhigh/fs;   % resample rate to get to receiver sample rate

% Generate the signal emitted by the target
[x, noise_bw] = generate_signal(Nsym, fsym, sps_high, sps, span, beta, ...
    show_plots); 

% Add proper delays that correspond to target and emitter locations
[y1, true_delays, true_tdoas, ranges] = add_delay(x, targetPos, refPos, ...
    fhigh, show_plots);

% Downsample the signals so they correspond to the receiver sample rate
% This will cause fractional delay offsets that will exist when correlating
y2 = lower_samp_rate(y1, resample_rate, show_plots);
targetPos
true_tdoas
avg_coords = [0;0];
MSE_coords = [0 0;0 0];
unique_avg = 0;
for nn = 1:Ntrials
    % Add noise at the proper SNR levels for free space path losses
    y3 = add_noise(y2, tx_pwr_dbm, noise_bw, fc, ranges, show_plots);

    % Estimate the delay using the received signals
    [tdoas, corr_mag_sq, peak_idxs, lags] = ...
        get_tdoa(y3, wlen, nstds, fs, show_plots);

    % Refine the tdoas using a super resolution algorithm
    sample_idxs = [peak_idxs-1;peak_idxs;peak_idxs+1];
    lags = [lags-1;lags;lags+1];
    numpairs = size(corr_mag_sq,2); % number of correlated pairs of receivers
    for ii = 1:numpairs
        corr_peak_samples(:,ii) = corr_mag_sq(sample_idxs(:,ii), ii);
    end
    [tdoas2, diffs] = refine_tdoa(lags, corr_peak_samples, fs, show_plots);
    % abs(true_tdoas-tdoas)./true_tdoas
    % abs(true_tdoas-tdoas2)./true_tdoas
%     tdoas2 = tdoas;

    % Feed the refined TDOAs to a localization algorithm
%     [coords, unique] = geo_lsq(refPos, tdoas2);
    [coords, unique] = geo_sphere_int(refPos, tdoas2);
    
    % Compute statistical performance metrics
    avg_coords = avg_coords + coords; 
    MSE_coords = MSE_coords + (targetPos - coords)*(targetPos-coords)';
    unique_avg = unique_avg + unique;
end

avg_coords = avg_coords/Ntrials;
bias_coords = avg_coords - targetPos;
MSE_coords = MSE_coords/Ntrials;
covar_coords = MSE_coords - bias_coords*bias_coords';
mse_coords = trace(MSE_coords);
unique_avg = unique_avg/Ntrials;
unique = unique_avg;
tdoas = tdoas2;
end


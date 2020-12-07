function [results] = get_single_emitter3(sim_params, transmitter_params, ...
        receiver_params, channel_params, loc_alg_params)

c = 299792458;  % speed of light m/s    
    
% extract data and assign to local variables
show_plots = sim_params.show_plots;
apply_calibration = sim_params.apply_calibration;
calib_path = sim_params.calib_path;
Ntrials = sim_params.Ntrials;

Nsym = transmitter_params.Nsyms;
fsym = transmitter_params.symbol_rate;
sps = transmitter_params.sps;
mod_type = transmitter_params.mod_type;
mod_order = transmitter_params.mod_order;
span = transmitter_params.span;
beta = transmitter_params.excess_bw;
fc = transmitter_params.carrier_freq;
tx_pwr_dbm = transmitter_params.tx_pwr_dbm;
targetPos = transmitter_params.target_locs;

fs = receiver_params.sample_rate;
refPos = receiver_params.ref_locs;
[numdims, numrefs] = size(refPos);
numpairs = numrefs - 1;         % unique pairs of receivers

multi_option = channel_params.multi_option;
delay_spread = channel_params.delay_spread;
num_paths = channel_params.num_paths;
max_num_paths = channel_params.max_num_paths;
multi_idxs = channel_params.multi_idxs;
max_nlos_amp = channel_params.max_amp;
min_num_taps = channel_params.min_num_taps;
multi_dist_based = channel_params.distance_based;
multi_coords = channel_params.multi_coords;
multi_delays = repmat({nan}, numrefs, 1);
if multi_option == 4
    for nn = 1:length(multi_coords)
        if ~isnan(multi_coords{nn})
            coords = multi_coords{nn};
            multi_delays{nn} = ( vecnorm(coords-refPos(:,nn)) + ...
                    vecnorm(coords-targetPos) )/c;
        end
    end
end

% grid = 0; % if dpd is not used it stays this by default

% This switch assigns struct variables to local variables depending on the
% localization algorithm used. This is separate to the section that runs
% the simulation because doing it this way minimizes code copying in that
% section
switch loc_alg_params.type
    case {'dpd','direct-dpd','ms-dpd1','ms-dpd2','ms-music-dpd','ms-music-dpd2','sparse-dpd'}
        grid_def = loc_alg_params.grid_def;
        
    case {'lsq', 'si'}
        percent_of_peak = receiver_params.percent_of_peak;
        unique_avg = 0;
        total_missed_peaks = 0;
        detection_count = 0;
        Ntrials_with_peak = Ntrials;    % if a peak is missed we skip a trial
        coords = nan(numdims, Ntrials);
        tdoas_coarse = nan(Ntrials, numpairs);
        tdoas_refined = nan(Ntrials, numpairs);
        tdoas_f = nan(Ntrials, numpairs);
        
    case 'taylor'
        initial_coords = loc_alg_params.initial_coords;
        percent_of_peak = receiver_params.percent_of_peak;
        unique_avg = 0;
        total_missed_peaks = 0;
        detection_count = 0;
        Ntrials_with_peak = Ntrials;    % if a peak is missed we skip a trial
        coords = nan(numdims, Ntrials);
        tdoas_coarse = nan(Ntrials, numpairs);
        tdoas_refined = nan(Ntrials, numpairs);
        tdoas_f = nan(Ntrials, numpairs);
end

% Generate the signal emitted by the target
% [x, noise_bw] = generate_signal2(Nsym, fsym, sps, span, beta, show_plots);
[x, noise_bw] = signal_generator(mod_type, Nsym, mod_order, sps, fc, ...
    fsym*sps, beta, span);
results.tx_raw = x;

[P,Q] = rat(fs/(fsym*sps));
y1 = resample(x, P, Q);
y1 = y1/sqrt(mean(abs(y1).^2)); % renormalize power to 1 after resample

% Set the transmitted signal amplitude to correspond to transmit signal power
tx_pwr = 10^(tx_pwr_dbm/10);
tx_amp = sqrt(tx_pwr);
y1 = tx_amp * y1;
results.tx_up = y1;

% Add proper delays that correspond to target and emitter locations
[y2, toas_true, tdoas_true, ranges] = add_delay2(y1, targetPos, refPos, ...
    fs, show_plots);
results.tx_delayed = y2;
results.ranges = ranges;
% tdoas_true

% Add proper attenuation for the LOS path signals
[y2p5] = add_attenuation(y2, fc, ranges);

avg_coords = [0;0];
MSE_coords = [0 0;0 0];
corr_mag_bad_sv = cell(1,numpairs);
kk = ones(1,numpairs);
for nn = 1:Ntrials
    
    % Add multipath
    [y3, h{nn}] = add_multipath(y2p5, y1, fc, fs, ranges, delay_spread, num_paths, ...
    max_num_paths, multi_idxs, multi_option, max_nlos_amp, min_num_taps, ...
    multi_dist_based, multi_delays, tx_pwr_dbm, show_plots);
    results.tx_multi{nn} = y3;
    
    % Add noise at the proper SNR levels for free space path losses
    [y4, avg_snr_db(:,nn)] = add_noise(y3, tx_pwr_dbm, noise_bw, fc, ranges, show_plots);
%     y4 = y3;
    results.tx_noise{nn} = y4;
    
    % Perform an AGC operation to received samples
    y5 = y4./max(abs(y4));
    results.tx_agc{nn} = y5;
    
%     figure
%     s1 = y1;
%     s2 = y4(:,1);
%     subplot(2,1,1)
%     plot(real(s1), 'bx-'); hold on
%     plot(real(s2), 'r.-')
%     subplot(2,1,2)
%     plot(imag(s1), 'bx-'); hold on
%     plot(imag(s2), 'r.-')

    switch loc_alg_params.type
        case {'dpd','direct-dpd','ms-dpd1','ms-dpd2','ms-music-dpd','ms-music-dpd2','sparse-dpd'}
            switch loc_alg_params.type 
                case 'dpd'
                    [coords(:,nn), grid, dpd_obj{nn}, unique] = dpd(y5, fs, refPos, grid_def);

                case 'direct-dpd'

        %             p0 = targetPos + 3*randn(2,1); % help the minimizer for testing
                    p0 = targetPos; % help the minimizer for testing
        %             p0 = [50;17.5];
                    [coords(:,nn), grid, dpd_obj{nn}, unique] = ...
                        dpd_direct(y4, fs, fc, refPos, grid_def, y1, tx_pwr_dbm, p0);

                case 'ms-dpd1'

                    [coords(:,nn), rcoords(:,nn), grid, reflect_grid, dpd_obj{nn}, unique] = ...
                        ms_dpd1(y4, fs, fc, refPos, grid_def, y1, tx_pwr_dbm);
                    
                case 'ms-dpd2'

                    [coords(:,nn), rcoords(:,nn), grid, reflect_grid, dpd_obj{nn}, unique] = ...
                        ms_dpd2(y4, fs, fc, refPos, grid_def, y1, tx_pwr_dbm);
                    
                case 'ms-music-dpd'
                    
                    [coords(:,nn), rcoords(:,nn), grid, reflect_grid, dpd_obj{nn}, unique] = ...
                        ms_music_dpd(y4, fs, fc, refPos, grid_def, y1, tx_pwr_dbm);
                    
                case 'ms-music-dpd2'

                [coords(:,nn), rcoords(:,nn), grid, reflect_grid, dpd_obj{nn}, unique] = ...
                    ms_music_dpd2(y4, fs, refPos, grid_def);
                
                case 'sparse-dpd'

                [coords(:,nn), all_coords{nn}, grid, sparse_heatmap{nn}, unique] = ...
                    sparse_dpd(y4, y1, fs, refPos, grid_def);
            end

            % Compute statistical performance metrics
            avg_coords = avg_coords + coords(:,nn); 
            MSE_coords = MSE_coords + (targetPos - coords(:,nn))*(targetPos-coords(:,nn))';

        case {'lsq', 'si', 'taylor'}
            % Estimate the delay using the received signals
            [tdoas_coarse(nn,:), tdoas_f(nn,:), corr_mag, peak_idxs, lags, ...
                lags_full, num_samps_from_peak, corr_mag_bad] = ...
                get_tdoa(y5, fs, percent_of_peak, apply_calibration, ...
                calib_path, tdoas_true, show_plots);
%             rx_num = 1;
            corr_mag_sv(nn,:,:) = corr_mag;
            for jj = 1:numpairs
                if ~isnan(corr_mag_bad{jj})
                    corr_mag_bad_sv{kk(jj),jj} = corr_mag_bad{jj};
                    kk(jj) = kk(jj) + 1;
                end
            end
            
            % for correlation plot debug
%             if nn == Ntrials
%             figure
%             rx_pair = 1;
%             plot(lags_full(:,1)/(fs*1e-9),squeeze(corr_mag_sv(:,:,rx_pair)));
%             end

            if sum(isnan(peak_idxs)) ~= numpairs
                detection_count = detection_count + 1;
            end
            peak_count = sum(~isnan(peak_idxs));
            missed_peak_count = sum(isnan(peak_idxs)); 
            total_missed_peaks = total_missed_peaks + missed_peak_count;

            if missed_peak_count == 0
                % Refine the tdoas using a super resolution algorithm
                numpairs = size(corr_mag,2); % number of correlated pairs of receivers  
                lags_new = zeros(2*num_samps_from_peak+1, numpairs);
                corr_peak_samples = zeros(2*num_samps_from_peak+1,numpairs);
                for ii = 1:numpairs
                    sample_idxs = (peak_idxs(ii)-num_samps_from_peak:peak_idxs(ii)+num_samps_from_peak).';
                    lags_new(:,ii) = (lags(ii)-num_samps_from_peak:lags(ii)+num_samps_from_peak).';
                    if isnan(sample_idxs)
                        corr_peak_samples(:,ii) = NaN;
                    else
                        corr_peak_samples(:,ii) = corr_mag(sample_idxs, ii);
                    end
                end
                corr_peak_samples = {corr_peak_samples; corr_mag};    % temporary while debugging sinc interp
                lags_new = {lags_new; lags_full};
                [tdoas_refined(nn,:), diffs] = refine_tdoa(lags_new, corr_peak_samples, fs, show_plots);

                % Feed the refined TDOAs to a localization algorithm
                switch loc_alg_params.type
                    case 'lsq'
                        [coords(:,nn), unique] = geo_lsq(refPos, tdoas_refined(nn,:));
        %               [coords(:,nn), unique] = geo_lsq(refPos, tdoas_f(nn,:));
                        
                    case 'si'
                        [coords(:,nn), unique] = geo_sphere_int(refPos, tdoas_refined(nn,:));
                        
                    case 'taylor'
                        [coords(:,nn), unique] = taylor_linearization(refPos, tdoas_refined(nn,:),initial_coords,[]);
                end                

                % Compute statistical performance metrics
                avg_coords = avg_coords + coords(:,nn); 
                MSE_coords = MSE_coords + (targetPos - coords(:,nn))*(targetPos-coords(:,nn))';
                unique_avg = unique_avg + unique;

            else
                Ntrials_with_peak = Ntrials_with_peak - 1;
                tdoas_refined(nn,:) = nan*ones(1,numpairs);
            end

        otherwise
            warning('Unexpected location algorithm type. Simulation skipped.')
    end
end

switch loc_alg_params.type
    case {'dpd','direct-dpd','ms-dpd1','ms-dpd2','ms-music-dpd','ms-music-dpd2','sparse-dpd'}
        if contains(loc_alg_params.type, 'ms-dpd')
            results.reflection_coords = rcoords;
            results.reflection_grid = reflect_grid;          
        end
        results.grid = grid;
        if ~contains(loc_alg_params.type, 'sparse-dpd')
            results.obj_vals = dpd_obj;
        end
        if contains(loc_alg_params.type, 'sparse-dpd')
            results.sparse_heatmap = sparse_heatmap;
            results.max_coord = coords;
            results.all_coords = all_coords;
        end
        avg_coords = avg_coords/Ntrials;
        MSE_coords = MSE_coords/Ntrials;

    case {'lsq', 'si', 'taylor'}
        avg_coords = avg_coords/Ntrials_with_peak;
        MSE_coords = MSE_coords/Ntrials_with_peak;
        unique_avg = unique_avg/Ntrials_with_peak;
        unique = unique_avg;
        prob_correlation = 1 - total_missed_peaks/(numpairs*Ntrials);
        prob_detection = detection_count/Ntrials;
        
        results.tdoas_coarse = tdoas_coarse;
        results.tdoas_refined = tdoas_refined;
        results.prob_correlation = prob_correlation;
        results.prob_detection = prob_detection;
        results.corr_mag = corr_mag_sv;
        results.corr_mag_bad = corr_mag_bad_sv;
        results.unique = unique;
end

bias_coords = avg_coords - targetPos;
covar_coords = MSE_coords - bias_coords*bias_coords';
mse_coords = trace(MSE_coords);

results.tdoas_true = tdoas_true;
results.toas_true = toas_true;
results.coords = coords;
results.bias_coords = bias_coords;
results.covar_coords = covar_coords;
results.mse_coords = mse_coords;
results.avg_snr_db = avg_snr_db;
results.channel_taps = h;
results.multi_delays = multi_delays;

end


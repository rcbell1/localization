function [coords, tdoas] = get_single_emitter_fromfile(file_path, ...
    refPos, wlen, nstds, percent_of_peak, show_plots)

[numdims, numrefs] = size(refPos);
numpairs = numrefs - 1;         % unique pairs of receivers

% Generate the signal emitted by the target
[y, fs] = load_signals_fromfile(file_path, show_plots);

coords = [inf;inf];

% Estimate the delay using the received signals
[tdoas, corr_mag_sq, peak_idxs, lags, num_samps_from_peak] = ...
    get_tdoa(y, wlen, nstds, fs, percent_of_peak, show_plots);

missed_peak_count = sum(isnan(peak_idxs)); 

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
    [tdoas_refined, diffs] = refine_tdoa(lags_new, corr_peak_samples, fs, show_plots);

    % Feed the refined TDOAs to a localization algorithm
    [coords, unique] = geo_lsq(refPos, tdoas_refined);
%         [coords, unique] = geo_sphere_int(refPos, tdoas2);
else
    tdoas_refined = nan*ones(1,numpairs);
end

tdoas = tdoas_refined;
end


function [tdoas_refined, tdoas_diff] = refine_tdoa(lags, corr_mag_sq, fs, show_plots)

% just for debug of sinc interp
corr_mag_full = corr_mag_sq{2};
corr_mag_sq = corr_mag_sq{1};
lags_full = lags{2};
lags = lags{1};

[npoints, numrefs] = size(lags);
coarse_peak_idx = (npoints+1)/2;    % center idx of mag_corr_sq array
Ts = 1/fs;
f = @(x,t) t(1)*x.^2 + t(2)*x + t(3);   % quadratic fit

if show_plots == 1
    hf1 = figure;
    hf2 = figure;
end

for ii = 1:numrefs
    w = lags(:,ii);
    z = corr_mag_sq(:,ii);
    
    % Parabolic Interpolation
    H = [w.^2 w ones(length(w),1)];
    theta = H\z; % the coefficients of the line fitting the data points

    wm = -theta(2)/(2*theta(1));        % value of x where y is max
    wm_t = wm*Ts;
    tdoas_diff(ii) = abs(w(2)*Ts-wm_t);
    tdoas_refined(ii) = wm_t;
    
    % Sinc Interpolation
    up_rate = 4;
    corr_up1 = resample(corr_mag_full(:,ii), up_rate, 1);
    [mval, midx] = max(corr_up1);    
    
    if show_plots == 1
        nrows = ceil(numrefs/2);
        figure(hf1)
        subplot(nrows, 2, ii);
        zm = f(wm,theta);                 % y max value
        xu = linspace(min(w), max(w), 100);
        xtu = xu*Ts*1e9;
        plot(w*Ts*1e9,z,'.-','markersize',14); hold all
        plot(xtu,f(xu,theta),'--'); % plot parabola
        plot(wm_t*1e9,zm,'o');  % mark max parapola
        plot(lags(coarse_peak_idx,ii)*Ts*1e9, corr_mag_sq(coarse_peak_idx,ii), 'bo','markersize',8)
        axis([min(w*Ts*1e9) max(w*Ts*1e9) 0.95*min(z) 1.05*max(z)])
        title(sprintf('Refined TDOA Ref %i and 1', ii+1))
        xlabel('TDOA (ns)')
        ylabel('|Correlator Out|^2')
        
        figure(hf2)
        subplot(nrows, 2, ii);
        lagsup = lags_full(1,ii):1/up_rate:lags_full(end,ii)+(up_rate-1)/up_rate;
        plot(lagsup.'*Ts*1e9, corr_up1, '--','markersize',14); hold all
        plot(lagsup(1:up_rate:end).'*Ts*1e9, corr_mag_full(:,ii), '.-','markersize',14)
        plot(lagsup(midx).'*Ts*1e9, mval, 'o')
        plot(lags(coarse_peak_idx,ii)*Ts*1e9, corr_mag_sq(coarse_peak_idx,ii), 'bo','markersize',8)
    end
end
    
end


function [tdoas_refined, tdoas_diff] = refine_tdoa(lags, corr_mag_sq, fs, show_plots)

[npoints, numrefs] = size(lags);
coarse_peak_idx = (npoints+1)/2;    % center idx of mag_corr_sq array
Ts = 1/fs;
f = @(x,t) t(1)*x.^2 + t(2)*x + t(3);   % quadratic fit

if show_plots == 1
    figure
end

for ii = 1:numrefs
    w = lags(:,ii);
    z = corr_mag_sq(:,ii);
    
    H = [w.^2 w ones(length(w),1)];
    theta = H\z; % the coefficients of the line fitting the data points

    wm = -theta(2)/(2*theta(1));        % value of x where y is max
    wm_t = wm*Ts;
    tdoas_diff(ii) = abs(w(2)*Ts-wm_t);
    tdoas_refined(ii) = wm_t;
    
    if show_plots == 1
        nrows = ceil(numrefs/2);
        subplot(nrows, 2, ii)
        zm = f(wm,theta);                 % y max value
        xu = linspace(min(w), max(w), 100);
        xtu = xu*Ts;
        plot(w*Ts,z,'.-','markersize',14); hold all
        plot(xtu,f(xu,theta),'--');
        plot(wm_t,zm,'o');
        plot(lags(coarse_peak_idx,ii)*Ts, corr_mag_sq(coarse_peak_idx,ii), 'bo','markersize',8)
        axis([min(w*Ts) max(w*Ts) 0.95*min(z) 1.05*max(z)])
        title(sprintf('Refined TDOA Ref %i and 1', ii+1))
        xlabel('TDOA (s)')
        ylabel('|Correlator Out|^2')
    end
end
    
end


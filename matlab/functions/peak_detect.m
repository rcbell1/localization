function [idxs, xmean, xstd] = peak_detect(x, wlen, nstds, show_plots)
% this is for positive valued series only, as is the case after a magnitude
% squared operation. wlen should be odd

xm = movmax(x,wlen);    % smooths the series to reject false peaks
[~,idx_max] = max(xm);  % index where the max inside window begins
idx_max = idx_max + floor(wlen/2);  % the center of the window is true max

% [val, idxs] = max(x); % simple but will result in false peaks when a true
                        % peak does not exist

% check that the peak is truely a global outlier peak and not local
xmean = mean(x);
xstd = std(x);

for ii = 1:size(x,2)
    if (x(idx_max(ii),ii) - xmean(ii)) > nstds*xstd(ii)
        idxs(ii) = idx_max(ii);
    else
        idxs(ii) = [];
    end
end

if show_plots == 1
    figure
    plot(x,'.-'); hold all
    xmin = inf; xmax = 0;
    for ii = 1:size(x,2)
        plot([1 length(x)], [nstds*xstd(ii) nstds*xstd(ii)], '--');
        plot(idxs(ii), x(idxs(ii),ii), 'o')
        if xmin > min(x(idxs(ii)-1,ii),x(idxs(ii)+1,ii))
            xmin = min(x(idxs(ii)-1,ii),x(idxs(ii)+1,ii));
        end
        if xmax < x(idxs(ii),ii)
            xmax = x(idxs(ii),ii);
        end
    end
    title('Cross Correlation Peak Location Identification')
    xlabel('Sample Number')
    ylabel('|Correlator Out|^2')
    axis([-inf inf 0 1.1*xmax])
    
    axes('position',[0.6 0.5 0.27 0.3])
    plot(x,'.-'); hold all
    for ii = 1:size(x,2)
        plot(idxs(ii), x(idxs(ii),ii), 'o')
    end
    axis([idxs(ii)-8 idxs(ii)+15 0.95*xmin 1.05*xmax])
    title('Zoom Peaks')
    xlabel('Sample Number')
    ylabel('|Correlator Out|^2')
end

end


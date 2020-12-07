close all; clear

% path = {'14 - wired 02_20_2020 1021/tx_center/rfs7/1/';
%         '14 - wired 02_20_2020 1021/tx_center/rfs7/2/';
%         '14 - wired 02_20_2020 1021/tx_center/rfs7/3/'};
 
% path = {'outdoors/3-03_07_2020/tx_center_sync/rfs100/1/'};
% path = {'indoors/1-03_11_2020/tx_center_sync/rfs100/1/'};
% path = {'desk/18/tx_center_sync/rfs100/10/'};
% path = {'desk/19/tx_center/rfs100/1/'};
% path = {'5 - wireless tests desk/tx_center/rfs9/6/'};
path = {'6 - long wired tests/tx_center/rfs9/6/'};
% path = {'4 - short wired tests/tx_center/rfs9/1/'};
% path = {'15/tx_center/rfs4/1/';
%         '15/tx_center/rfs4/2/';
%         '15/tx_center/rfs4/3/'};
file_name = 'rx_pulses_wired.mat';
beta = 0.5;

num_files = length(path);
for jj = 1:num_files
    load([path{jj} file_name])
    Nsamp = double(Nprx);
    fs = fs_rx;
    Ntrials = length(bounds)-1;
    Nrx = size(yblock,2);

    %% Process all pulses in a file
    Nfft = 2^(nextpow2(2*Nsamp-1));
    faxis = fs*(-0.5:1/Nfft:0.5-1/Nfft);
    bw = fs/sps_rx*(1+beta/3); % the bandwidth to use for computations
    array_idx1 = find(faxis > -bw/2 & faxis < bw/2);
    removal_band = 0.07e6; % band around zero to ignore, no delay information at k=0
%     removal_band = 1e6;
    array_idx2 = find(faxis > -removal_band & faxis < removal_band);
    [~, shared_idx] = intersect(array_idx1,array_idx2,'stable');
    array_idx1(shared_idx) = []; % delete those frequency indices in the removal band
    array_idx = array_idx1;
    kidx = faxis*Nfft/fs;
    fkeep_idx = find(faxis > -bw/2 & faxis < bw/2);
    frem_idx = find(faxis > -removal_band & faxis < removal_band);
    [~, shared_f] = intersect(fkeep_idx,frem_idx,'stable');
    fkeep_idx(shared_idx) = []; % delete those frequency indices in the removal band
    fk = faxis(fkeep_idx);
    kidx = kidx(fkeep_idx);
    H = [fk.' ones(length(fk),1)];
    zero_idx = ceil(Nfft/2) + 1; % phase array index of frequency zero
        
    phases_f = zeros(Ntrials,Nfft, Nrx-1);
    phases_f_unwrap = zeros(Ntrials,Nfft,Nrx-1);
    corrs_t = zeros(Ntrials, 2*Nsamp-1,Nrx-1);
    corrs_f = zeros(Ntrials,Nfft, Nrx-1);
    for nn = 1:Ntrials
        % get frequency domain correlation and phases for each pulse
        y = yblock(bounds(nn):bounds(nn+1)-1,:);
        y_f = fft(y,Nfft);
%         Nseg = 2^13;
%         pxx = pwelch(y, Nseg, [], Nfft);
%         y_f = pxx;
        corrs_f(nn,:,:) = y_f(:,2:Nrx).*conj(y_f(:,1));
        
        % get the phase in frequency domain for each stream
        phases_f(nn,:,:) = angle(fftshift(corrs_f(nn,:,:),2));
        phases_f(nn,:,:) = squeeze(phases_f(nn,:,:)) - squeeze(mean(phases_f(nn,array_idx,:),2)).';
%         phases_f_unwrap(nn,array_idx,:) = unwrap(phases_f(nn,array_idx,:),[],2);
        phases_f_unwrap(nn,:,:) = phases_f(nn,:,:);
        phases_f_unwrap(nn,array_idx,:) = unwrap(phases_f(nn,array_idx,:),pi,2);
        phases_f_unwrap(nn,:,:) = squeeze(phases_f_unwrap(nn,:,:)) - squeeze(mean(phases_f_unwrap(nn,array_idx,:),2)).';
        
        % do a linear fit on the raw phase
        p = squeeze(phases_f(nn,array_idx,:));
        m = H\p;
        slopes = m(1,:);
        intercepts = m(2,:);
        p_est(nn,:,:) = repmat(fk.',1,Nrx-1)*diag(slopes) + intercepts;

        delays_est1(nn,:,:) = -Nfft*p./(2*pi*fs*kidx.')/1e-9;
        delays_est2(nn,:) = mean(-Nfft*squeeze(p_est(nn,:,:))./(2*pi*fs*kidx.'))/1e-9;

        % do a linear fit on the unwrapped phase        
        p_unwrap = squeeze(phases_f_unwrap(nn,array_idx,:));
        m = H\p_unwrap;
        slopes_unwrap = m(1,:);
        intercepts_unwrap = m(2,:);
        p_unwrap_est(nn,:,:) = repmat(fk.',1,Nrx-1)*diag(slopes_unwrap) + intercepts_unwrap;

        delays_unwrap_est1(nn,:,:) = -Nfft*p_unwrap./(2*pi*fs*kidx.')/1e-9;
        delays_unwrap_est2(nn,:) = mean(-Nfft*squeeze(p_unwrap_est(nn,:,:))./(2*pi*fs*kidx.'))/1e-9;
        
        % time domain correlation for each pulse
        for mm = 1:Nrx-1
            [corrs_t(nn,:,mm), lags] = xcorr(y(:,mm+1),y(:,1));
        end
    end
    mean_delays_est1(jj,:) = squeeze(mean(delays_est1,[1 2]))
    mean_delays_est2(jj,:) = mean(delays_est2)
    mean_delays_unwrap_est1(jj,:) = squeeze(mean(delays_unwrap_est1,[1 2]))
    mean_delays_unwrap_est2(jj,:) = mean(delays_unwrap_est2)
%     std_delays_est(jj,:) = std(delays_est);
    
    
    %% Plots for each file
    % Plot the time domain correlation between receivers
%     plot_idxs = [1 round(num_trials/2) num_trials];
    plot_idxs = 1;
    plot_opts = {'b-','r-','g-','o-'};
    
    figure
    subplot(2,3,[1,4])
    for nn = 1:length(plot_idxs)
        for mm = 1:Nrx-1
            plot(lags, abs(corrs_t(plot_idxs(nn),:,mm)), plot_opts{mm}, 'markersize', 10); hold on
        end
    end
    ylims = ylim;
    plot([0 0], [ylims(1) ylims(2)], 'k--')
    title('Time Domain Correlations')
    legend('C_{21}', 'C_{31}')
    xlabel('Sample Lag')
    axis([-10 10 -inf inf])

    % Plot the frequency domain phases for each receiver stream
    % Wrapped phases
    subplot(2,3,2)
    for nn = 1:length(plot_idxs)
        for mm = 1:Nrx-1
            yaxis = squeeze(phases_f(plot_idxs(nn),:,mm));
            f(mm) = plot(faxis/1e6, yaxis, plot_opts{mm}); hold on
%             plot(faxis(array_idx)/1e6, squeeze(phases_f(plot_idxs(nn),array_idx,mm)), 'k.')
            plot(faxis(array_idx)/1e6, p_est(plot_idxs(nn),:,mm), 'k--')
%             text(faxis(1)*0.95/1e6,-1.45*pi+(mm-1)*pi/7,sprintf('Estimated Delay: %3.1f ns', mean_delays_est2(Nrx-mm)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment','left','clipping','on','fontweight','bold')
        end
    end
    ylims = ylim;
    plot([-bw/2 -bw/2]/1e6, [-2*pi 2*ylims(2)], 'k--')
    plot([bw/2 bw/2]/1e6, [-2*pi 2*ylims(2)], 'k--')
    axis([-inf inf -1.5*pi ylims(2)])
%     title(sprintf('File %i - Spectral Phase', jj))
    title(sprintf('Spectral Phase - Indoor', jj))
    xlabel('Frequency (MHz)')
    ylabel('Phase')
    legend(f, 'C_{21}', 'C_{31}')

    % Unwrapped phases
    subplot(2,3,3)
    for nn = 1:length(plot_idxs)
        for mm = 1:Nrx-1
            yaxis = squeeze(phases_f_unwrap(plot_idxs(nn),:,mm));
            f(mm) = plot(faxis/1e6, yaxis, plot_opts{mm}); hold on
%             plot(faxis(array_idx)/1e6, squeeze(phases_f_unwrap(plot_idxs(nn),array_idx,mm)), 'k.')
            plot(faxis(array_idx)/1e6, p_unwrap_est(plot_idxs(nn),:,mm), 'k-')
        end
    end
    ylims = ylim;
    for mm = 1:Nrx-1
%         text(faxis(1)*0.95/1e6,ylims(2)-(mm-1)*.1*ylims(2),sprintf('Estimated Delay: %3.1f ns', mean_delays_unwrap_est2(mm)), 'VerticalAlignment', 'top', 'HorizontalAlignment','left','clipping','on','fontweight','bold')
    end
    plot([-bw/2 -bw/2]/1e6, [ylims(1) ylims(2)], 'k--')
    plot([bw/2 bw/2]/1e6, [ylims(1) ylims(2)], 'k--')
%     title('Spectral Phase (Unwrapped)')
    title('Spectral Phase - Outdoor')
    xlabel('Frequency (MHz)')
    ylabel('Phase')
    axis([-inf inf ylims(1) ylims(2)])
    legend(f, 'C_{21}', 'C_{31}')
   
    % Plot the phases converted to TDOA estimates
    % Wrapped phases
    plot_opts = {'b.','r.','g.','c.','y.','m.'};   
    subplot(2,3,5)
    for nn = 1:length(plot_idxs)
        for mm = 1:Nrx-1
            yaxis = squeeze(delays_est1(plot_idxs(nn),:,mm));
            f(mm) = plot(faxis(array_idx)/1e6, yaxis, plot_opts{mm}); hold on
            plot([-bw/2 bw/2]/1e6, [mean_delays_est1(mm) mean_delays_est1(mm)], 'k--')
            text((bw/2)/1e6,mean_delays_est1(mm),sprintf('Estimated Delay: %3.1f ns', mean_delays_est1(mm)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment','right','fontweight','bold')
        end
    end
    title('Per Bin Delay Estimate')
    xlabel('Frequency (MHz)')
    ylabel('Delay (ns)')
    legend(f, 'C_{21}', 'C_{31}')
    yb1 = min(mean_delays_est1);
    yb2 = max(3,max(mean_delays_est1));
    axis([-inf inf yb1-abs(yb1) 2*yb2])

    % Unwrapped phases
    subplot(2,3,6)
    for nn = 1:length(plot_idxs)
        for mm = 1:Nrx-1
            yaxis_u = squeeze(delays_unwrap_est1(plot_idxs(nn),:,mm));
            f(mm) = plot(faxis(array_idx)/1e6, yaxis_u, plot_opts{mm}); hold on
            plot([-bw/2 bw/2]/1e6, [mean_delays_unwrap_est1(mm) mean_delays_unwrap_est1(mm)], 'k--')
            text((bw/2)/1e6,mean_delays_unwrap_est1(mm),sprintf('Estimated Delay: %3.1f ns', mean_delays_unwrap_est1(mm)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment','right','fontweight','bold')
        end
    end
    title('Per Bin Delay Estimate (Unwrapped)')
    xlabel('Frequency (MHz)')
    ylabel('Delay (ns)')
%     axis([-inf inf -xb xb])
    legend(f, 'C_{21}', 'C_{31}')
    yb1 = min(mean_delays_unwrap_est1);
    yb2 = max(3,max(mean_delays_unwrap_est1));
    axis([-inf inf yb1-abs(yb1) 2*yb2])

end
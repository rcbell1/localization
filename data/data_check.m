close all; clear

% path = {'14 - wired 02_20_2020 1021/tx_center/rfs7/1/';
%         '14 - wired 02_20_2020 1021/tx_center/rfs7/2/';
%         '14 - wired 02_20_2020 1021/tx_center/rfs7/3/'};
 

path = {'15/tx_center/rfs7/1/';
        '15/tx_center/rfs7/2/';
        '15/tx_center/rfs7/3/'};
file_name = 'rx_pulses_sliced.mat';
rolloff = 0.5;

num_files = length(path);
for jj = 1:num_files
    load([path{jj} file_name])
    num_trials = length(bounds)-1;
    nrx = size(yblock,2);

    %% Process all pulses in a file
    nfft = 2*Nprx-1;
    phases_f = zeros(num_trials,nfft-1, nrx);
    phases_f_diffs = zeros(num_trials,nfft-1,nrx-1);
    corrs_t = zeros(num_trials, nfft,nrx-1);
    for nn = 1:num_trials
        % get frequency domain for each pulse
        y = yblock(bounds(nn):bounds(nn+1)-1,:);
        y_f = fft(y,nfft);

        % time domain correlation for each pulse
        [corrs_t(nn,:,1)] = xcorr(y(:,2),y(:,1));
        [corrs_t(nn,:,2), lags] = xcorr(y(:,3),y(:,1));
        
        % equivalent of correlation in freq domain
    %     yyc_f = y_f(:,2:3).*conj(y_f(:,1));
    %     yyc_f_avg = yyc_f_avg + yyc_f;  % get the avg value across all pulses

        % get the phase in frequency domain for each stream
        phases_f(nn,:,:) = unwrap(angle(fftshift(y_f(2:end,:),1)));
        phases_f_diffs(nn,:,:) = phases_f(nn,:,2:3) - phases_f(nn,:,1);
    end

    %% Plots for each file
    % Plot the time domain correlation between receivers
%     plot_idxs = [1 round(num_trials/2) num_trials];
%     plot_idxs = 1:ceil(num_trials/2);
    plot_idxs = 1:num_trials;
    plot_opts = {'b.-','r.-','g.-','o.-'};
    
    figure
    subplot(2,2,[1,3])
    for nn = 1:length(plot_idxs)
        for mm = 1:nrx-1
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
    bw = fs_tx/sps*(1+rolloff);
    plot_opts = {'b-','r-','g-','o-'};
    fbin_idxs = fftshift(1:nfft-1);
    divider = 2*pi*fbin_idxs;
%     divider((nfft+1)/2) = 1;
    faxis1 = fs_rx*(-0.5:1/nfft:0.5-1/nfft); 
    faxis = fs_rx*(-0.5:1/nfft:0.5-1/nfft);
    faxis((nfft+1)/2) = [];
%     figure
    subplot(2,2,2)
    for nn = 1:length(plot_idxs)
        for mm = 1:nrx
            plot(faxis/1e6, squeeze(phases_f(plot_idxs(nn),:,mm)./divider), plot_opts{mm}); hold on
        end
    end
    ylims = ylim;
    plot([-bw/2 -bw/2]/1e6, [ylims(1) ylims(2)], 'k--')
    plot([bw/2 bw/2]/1e6, [ylims(1) ylims(2)], 'k--')
    title(sprintf('File %i - Phase of Frequency Domain Samples', jj))
    xlabel('Frequency (MHz)')
    ylabel('Phase')
    legend('Rx1','Rx2','Rx3')

    % Plot the phases differences between receivers 2-1 and 3-1
    subplot(2,2,4)
    for nn = 1:length(plot_idxs)
        for mm = 1:nrx-1
            plot(faxis/1e6, squeeze(phases_f_diffs(plot_idxs(nn),:,mm)), plot_opts{mm}); hold on
        end
    end
    ylims = ylim;
    plot([-bw/2 -bw/2]/1e6, [ylims(1) ylims(2)], 'k--')
    plot([bw/2 bw/2]/1e6, [ylims(1) ylims(2)], 'k--')
    title('Phase Differences of Frequency Domain Samples')
    xlabel('Frequency (MHz)')
    ylabel('Phase Difference')
    legend('Rx21','Rx31')
end
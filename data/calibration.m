clear; close all

file_path = '14 - wired 02_20_2020 1021/tx_center/rfs7/1/';
file_path = '4 - short wired tests/tx_center/rfs9/1/';
file_name = 'rx_pulses_sliced.mat';

load([file_path file_name])
fs = fs_rx;
num_trials = length(bounds)-1;
[nsamp_pulse, nrx] = size(yblock(bounds(1):bounds(2)-1,:));

%% For debug only, add false but known delays to input
% yblock = repmat(yblock(:,1),1,nrx);
% delays = [0 10e-9 250e-9]; % delays for each receiver (s)
% max_num_delay_samps = ceil(max(delays)*fs);
% for ii = 1:nrx
%     yblock_new(:,ii) = delayseq([yblock(:,ii); zeros(max_num_delay_samps,1)], delays(ii), fs);
% end
% yblock = yblock_new;
% end debug

%% Compute calibration values
nfft = 2*nsamp_pulse-1;
yyc_f_avg = zeros(nfft, nrx-1);
phases_f = zeros(num_trials,nfft, nrx);
phases_f_diffs = zeros(num_trials,nfft,nrx-1);
for nn = 1:num_trials
    % get frequency domain for each pulse
    y = yblock(bounds(nn):bounds(nn+1)-1,:);
    y_f = fft(y,nfft);
    
    % equivalent of correlation in freq domain
    yyc_f = y_f(:,2:3).*conj(y_f(:,1));
    yyc_f_avg = yyc_f_avg + yyc_f;  % get the avg value across all pulses
    
    % get the phase in frequency domain for each stream
    phases_f(nn,:,:) = unwrap(angle(y_f));
    phases_f_diffs(nn,:,:) = phases_f(nn,:,2:3) - phases_f(nn,:,1);
end
yyc_f_avg = yyc_f_avg/num_trials; % calibrate each pulse with this

save([file_path 'calibration.mat'], 'yyc_f_avg');

%% Debug, test the calibration estimate
% prove xcorr in freq matches xcorr in time
figure
lags_f = (1:nfft)-nsamp_pulse;
corr_f = fftshift(ifft(y_f(:,2).*conj(y_f(:,1))),1);
[corr_t, lags_t] = xcorr(y(:,2),y(:,1));
plot(lags_f, abs(corr_f), 'b.-', 'markersize', 10); hold on
plot(lags_t, abs(corr_t), 'ro-')
axis([-60 60 -inf inf])
title('Comparing Freq Domain Xcorr to Time Domain for Sanity')
xlabel('Sample Lag')
legend('Freq Domain','Time Domain')

% Show the xcorr of raw input pulse before calibration
figure
subplot(2,1,1)
[y_nc12, lags_nc12] = xcorr(y(:,2),y(:,1)); hold on
[y_nc13, lags_nc13] = xcorr(y(:,3),y(:,1));

plot(lags_nc12, abs(y_nc12), 'b.-', 'markersize', 10)
plot(lags_nc13, abs(y_nc13), 'r.-', 'markersize', 10)
title('No Calibration Applied')
legend('C_{21}', 'C_{31}')
xlabel('Sample Lag')
axis([-10 10 -inf inf])

subplot(2,1,2)
lags = (1:nfft)-nsamp_pulse;
y_c_avg = fftshift(ifft(yyc_f.*conj(yyc_f_avg)),1);
y_c_last = fftshift(ifft(yyc_f.*conj(yyc_f)),1);
plot(lags, abs(y_c_avg(:,1)), 'b.-', 'markersize', 10); hold on
plot(lags, abs(y_c_last(:,1)), 'ro-')
title(sprintf('Comparing estimated calibration of single pulse using \nall pulses to self calibration of single pulse with itself for sanity'))
legend('Estimated Avg Calibration', 'Self Signal Calibration')
xlabel('Sample Lag')
axis([-10 10 -inf inf])

% Plot the frequency domain phases for each of the receiver streams
faxis = fs_rx*(-0.5:1/nfft:0.5-1/nfft);
idxs = [1 round(num_trials/2) num_trials];
plot_opts = {'b-','r-','g-','o-'};
figure
subplot(2,1,1)
for nn = 1:length(idxs)
    for mm = 1:nrx
        plot(faxis/1e6, squeeze(phases_f(nn,:,mm)), plot_opts{mm}); hold on
    end
end
ylims = ylim;
plot([-fs_tx/2 -fs_tx/2]/1e6, [ylims(1) ylims(2)], 'k--')
plot([fs_tx/2 fs_tx/2]/1e6, [ylims(1) ylims(2)], 'k--')
title('Phase of Frequency Domain Samples')
xlabel('Frequency (MHz)')
ylabel('Phase')
legend('Rx1','Rx2','Rx3')

subplot(2,1,2)
for nn = 1:length(idxs)
    for mm = 1:nrx-1
        plot(faxis/1e6, squeeze(phases_f_diffs(nn,:,mm)), plot_opts{mm}); hold on
    end
end
ylims = ylim;
plot([-fs_tx/2 -fs_tx/2]/1e6, [ylims(1) ylims(2)], 'k--')
plot([fs_tx/2 fs_tx/2]/1e6, [ylims(1) ylims(2)], 'k--')
title('Phase Differences of Frequency Domain Samples')
xlabel('Frequency (MHz)')
ylabel('Phase Difference')
legend('Rx21','Rx31')
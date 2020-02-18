% This script demonstrates how you can use the fft and ifft to deconvolve
% sequences. 
%
% There are subutlies here between circular convolution and linear
% convolution that have to be understood. When you depend on the DFT to do
% frequency domain equalization, you are inherently assuming a circular
% convolution took place in the time domain because division of one DFT by 
% another only undoes circular convolution. If this is not true, your
% estimation using it will be wrong. There are two ways to enforce the
% circular convolution necessity. 
%   1) If you are constrained to use a predetermined length FFT/IFFT pair, 
%      as is the case when implementing OFDM systems, then you must 
%      implement the cyclic prefix concept into your signal structure. 
%   2) If you are not constrained by a fixed FFT/IFFT size, then you can 
%      make sure your FFT length is at least as long as N1+N2-1, where N1 
%      is the length of the first sequence and N2 the second. Without using
%      a cyclic prefix but maintaining the minimum FFT size you are
%      ensuring that linear convolution and circular convolution produce 
%      the same result.
clear; close all

path = '14/tx_center/rfs9/1/';
file_rx = 'rx_pulses_sliced.mat';
file_tx_pulse = 'tx_pulse.dat';
file_tx_pulse_clean = 'tx_pulse_clean.dat';
file_tx_sym = 'tx_symbols.dat';
load([path file_rx])
spreamble = read_complex_binary([path file_tx_pulse]);
spreamble_clean = read_complex_binary([path file_tx_pulse_clean]);
preamble = read_complex_binary([path file_tx_sym]);

[P,Q] = rat(fs_rx/fs_tx);
spreamble_rs = resample(spreamble_clean, P, Q);

N1 = length(spreamble); % length of preamble sequence in samples
delay_spread = 300;     % channel delay spread (ns)

tap_thresh = 0.01;  % magnitude cutoff for channel tap length estimate
span = 10;          % span in symbols of pulse shape filter
beta = 0.5;         % roll off factor of pulse shape filter
% sps = 4;            % samples per symbol
delta = 1e-8;       % channel tap error counting distance 

N2 = ceil(fs_rx*delay_spread*1e-9); % length of channel sequence in samples
Nout = 1*(N1 + N2 - 1); % the length of the channel output sequence

% preamble = 2*randi([0 1], N1, 1)-1 +1j*(2*randi([0 1], N1, 1)-1);
% rrc = rcosdesign(beta, span, sps);
% rrc = rrc/max(rrc);
% spreamble = conv(rrc, upsample(preamble,sps)).';
% % spreamble = spreamble(span+1:end-span); % strip garbage samples 
% N11 = length(spreamble); % length of shaped preamble in samples
% Nout = 1*(N11 + N2 - 1);  % length of linear convolution in samples
% % channel = [0.8 0 -0.1 0 0 0.1-1j*0.01];
% % channel = exponential(N2, 2);
% channel = [1 zeros(1,N2-1)];
% 
% received = conv(spreamble, channel);

% now suppose we are given the first sequence and the convolution output
% and we want to recover the second sequence
% use the known preamble to find start of received preamble

for nn = 1:length(bounds)-1
    y = yblock(bounds(nn):bounds(nn+1)-1,:);
    [corr1, lags1] = xcorr(y(:,1), spreamble_rs);
    mag_corr1 = abs(corr1);
    [mval, midx] = max(mag_corr1);
    Ndelay_est = lags1(midx);
    received(nn,:) = y(Ndelay_est+1:Ndelay_est+6*Nout,1);
    received(nn,:) = received(nn,:)/max(received(nn,:));
%     received(nn,:) = y(:,1);

    ffty = fft(received(nn,:), Nout);
    fftx = fft(spreamble, Nout);
    channel_est_f = ffty./fftx;
    channel_est = ifft(channel_est_f);
    num_taps_est = find(abs(channel_est)<tap_thresh,1)-1;
%     nerrors(nn) = sum(abs(channel - channel_est(1:N2)) > delta);
end
% avg_nerrors = mean(nerrors);
% nerrors = sum(abs(channel - channel_est(1:N2)) > delta);

% fprintf(1, '\n\nNumber of channel tap estimate errors on average: %i\n\n', avg_nerrors)

%% Plots
figure
plot(lags1, mag_corr1)

%% Plots
figure
subplot(4,2,1)
stem(real(preamble), 'filled')
xlabel('Symbol Number')
ylabel('Amplitude')
title('Real Preamble Symbols')
axis([-inf inf -1.5 1.5])
grid on

subplot(4,2,2)
stem(imag(preamble), 'filled')
xlabel('Symbol Number')
ylabel('Amplitude')
title('Imaginary Preamble Symbols')
axis([-inf inf -1.5 1.5])
grid on

subplot(4,2,3)
num_junk = span*sps/2+1; 
plot(1:N1, real(spreamble), '.-'); hold all
plot(num_junk:sps:N1-num_junk, real(preamble).', '.', 'markersize', 12)
xlabel('Sample Number')
ylabel('Amplitude')
title('Real Shaped Transmitted Preamble')
axis([1 400 -1.5 1.5])
grid on

subplot(4,2,4)
plot(1:N1, imag(spreamble), '.-'); hold all
plot(num_junk:sps:N1-num_junk, imag(preamble), '.', 'markersize', 12)
xlabel('Sample Number')
ylabel('Amplitude')
title('Imaginary Shaped Transmitted Preamble')
axis([1 400 -1.5 1.5])
grid on

% subplot(5,2,5)
% stem(real(channel), 'filled')
% xlabel('Tap Number')
% ylabel('Amplitude')
% title('Real Channel Response')
% axis([-inf inf -1.1 1.1])
% grid on
% 
% subplot(5,2,6)
% stem(imag(channel), 'filled')
% xlabel('Tap Number')
% ylabel('Amplitude')
% title('Imaginary Channel Response')
% axis([-inf inf -1.1 1.1])
% grid on

subplot(4,2,5)
num_junk = round(span*sps_rx/2+1); 
xaxis = num_junk:round(sps_rx):N1-num_junk;
xaxis = 1:round(sps_rx):length(spreamble);
len_xaxis = length(xaxis);
plot(real(received(end-1,:)), '.-'); hold all
plot(xaxis, real(preamble(1:len_xaxis)), '.', 'markersize', 12)
xlabel('Sample Number')
ylabel('Amplitude')
title('Real Received Preamble')
axis([1 400 -1.5 1.5])
axis([1 inf -1.5 1.5])
grid on

subplot(4,2,6)
xaxis = 1:round(sps_rx):length(spreamble);
plot(imag(received(end-1,:)), '.-'); hold all
plot(xaxis, imag(preamble(1:len_xaxis)), '.', 'markersize', 12)
xlabel('Sample Number')
ylabel('Amplitude')
title('Imaginary Received Preamble')
axis([1 400 -1.5 1.5])
axis([1 inf -1.5 1.5])
grid on

subplot(4,2,7)
stem(real(channel_est(1:num_taps_est)), 'filled')
xlabel('Tap Number')
ylabel('Amplitude')
title('Real Channel Estimate')
axis([-inf inf -5 5])
grid on

subplot(4,2,8)
stem(imag(channel_est(1:num_taps_est)), 'filled')
xlabel('Tap Number')
ylabel('Amplitude')
title('Imaginary Channel Estimate')
axis([-inf inf -5 5])
grid on

%% Plots
% For plotting purposes, the Nfft length needs to be a multiple of
% channel_est_f length, otherwise vector length problems will occur.
figure
fft_up = 256;
Nfft = length(channel_est_f)*fft_up;
faxis = -0.5:1/Nfft:0.5-1/Nfft;
% channel_f = abs(fftshift(fft(channel,Nfft)));
% plot(faxis, channel_f); hold on

% For some reason if N2, the channel tap length, is even, I have to start
% at the first sample after half the upsample rate of the fft to get the
% correct points. i.e. if N2 is even, then the first faxis index should be
% fft_up/2+1. If N2 is odd, then the first faxis index should be 1.
if mod(N2,2) == 0
    plot(faxis(129:fft_up:end), abs(fftshift(channel_est_f)).', 'r.')
else
    plot(faxis(1:fft_up:end), abs(fftshift(channel_est_f)).', 'r.')
end
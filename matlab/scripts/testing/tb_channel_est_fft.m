% This script demonstrates how you can use the fft and ifft to deconvolve
% sequences. 
%
% There are subtlties here between circular convolution and linear
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

% delay_spread = 300;     % channel delay spread (ns)

N1 = 10;            % length of preamble
span = 10;          % span in symbols of pulse shape filter
beta = 0.4;         % roll off factor of pulse shape filter
sps = 2;            % samples per symbol
% delta = 1e-8;       % channel tap error counting distance 

% N2 = ceil(fs_rx*delay_spread*1e-9); % length of channel sequence in samples
% Nout = 1*(N1 + N2 - 1); % the length of the channel output sequence


preamble = 2*randi([0 1], N1, 1)-1 + 1j*(2*randi([0 1], N1, 1)-1);
rrc = rcosdesign(beta, span, sps);
rrc = rrc/max(rrc);
spreamble = conv(rrc, upsample(preamble,sps)).'; % shaped preamble
% spreamble = spreamble(span+1:end-span); % strip garbage samples 
N1s = length(spreamble); % length of shaped preamble in samples
Nout = 1*(2*N1s);  % length of linear convolution in samples
channel = [zeros(1,5) 0.8 0 -0.3+0.4j 0 0 0.5-1j*0.3];
% channel = exponential(N2, 2);
% channel = 1;

received = conv(spreamble, channel);

channel_est = channel_est_fft(spreamble, received);

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
plot(1:N1s, real(spreamble), '.-'); hold all
plot(num_junk:sps:N1s-num_junk, real(preamble).', '.', 'markersize', 12)
xlabel('Sample Number')
ylabel('Amplitude')
title('Real Shaped Transmitted Preamble')
axis([1 inf -1.5 1.5])
grid on

subplot(4,2,4)
plot(1:N1s, imag(spreamble), '.-'); hold all
plot(num_junk:sps:N1s-num_junk, imag(preamble), '.', 'markersize', 12)
xlabel('Sample Number')
ylabel('Amplitude')
title('Imaginary Shaped Transmitted Preamble')
axis([1 inf -1.5 1.5])
grid on

subplot(4,2,5)
% num_junk = round(span*sps/2+1); 
% xaxis = num_junk:round(sps):N1s-num_junk;
xaxis = 1:round(sps):length(spreamble);
xaxis = 1:length(preamble);
len_xaxis = length(xaxis);
plot(real(received), '.-'); hold all
% plot(xaxis, real(preamble), '.', 'markersize', 12)
xlabel('Sample Number')
ylabel('Amplitude')
title('Real Received Preamble')
axis([1 inf -1.5 1.5])
axis([1 inf -1.5 1.5])
grid on

subplot(4,2,6)
xaxis = 1:round(sps):length(spreamble);
plot(imag(received), '.-'); hold all
% plot(xaxis, imag(spreamble(1:len_xaxis)), '.', 'markersize', 12)
xlabel('Sample Number')
ylabel('Amplitude')
title('Imaginary Received Preamble')
axis([1 inf -1.5 1.5])
grid on

subplot(4,2,7)
stem(real(channel_est(1:end)), 'filled', 'MarkerSize', 3)
xlabel('Tap Number')
ylabel('Amplitude')
title('Real Channel Estimate')
ylims = [-1 1];
axis([-inf inf ylims])
grid on

subplot(4,2,8)
stem(imag(channel_est(1:end)), 'filled', 'MarkerSize', 3)
xlabel('Tap Number')
ylabel('Amplitude')
title('Imaginary Channel Estimate')
ylims = [-1 1];
axis([-inf inf ylims])
grid on

%% Plots
% For plotting purposes, the Nfft length needs to be a multiple of
% channel_est_f length, otherwise vector length problems will occur.
figure
% fft_up = 256;
Nfft = 1024;
channel_est_f = fft(channel_est, Nfft);
% Nfft = length(channel_est_f)*fft_up;
faxis = -0.5:1/Nfft:0.5-1/Nfft;
plot(faxis, abs(fftshift(channel_est_f)).', 'r.')
% channel_f = abs(fftshift(fft(channel,Nfft)));
% plot(faxis, channel_f); hold on

% For some reason if N2, the channel tap length, is even, I have to start
% at the first sample after half the upsample rate of the fft to get the
% correct points. i.e. if N2 is even, then the first faxis index should be
% fft_up/2+1. If N2 is odd, then the first faxis index should be 1.
% if mod(N2,2) == 0
%     plot(faxis(129:fft_up:end), abs(fftshift(channel_est_f)).', 'r.')
% else
%     plot(faxis(1:fft_up:end), abs(fftshift(channel_est_f)).', 'r.')
% end
clear; close all

N = 30;     % length of first sequence (symbols)
fs = 10e6;  % receiver sampling rate (Hz)
fbw = 7e6;  % occupied bandwidth of transmitter
Ts = 1/fs;  % sample period (s)
D = 3*(2*rand-1)*Ts;  % delay (s)
sig = 0.05; % standard deviation of noise
sps1 = 5;  % samples per sample at high rate when multipath is added
sps2 = 4;   % samples per symbol for shaping filter
M = 4;      % modulation order
span = 5;   % length of shaping filter (symbols)
beta = 0.2; % excess bandwidth of shaping filter

% create the baseline modulated signal
xb = randi([0 M-1], log2(M)*N, 1);
x1 = qammod(xb, M, 'gray');
x2 = upsample(x1, sps2);
rrc = rcosdesign(beta, span, sps2, 'sqrt');
rrc = rrc.'/max(rrc);
x3 = conv(rrc, x2);
[P,Q] = rat((fs/fbw)/sps2);
x4 = resample(x3, P, Q);
% x4=x3;
x5 = x4/sqrt(mean(abs(x4).^2)); % normalize power to 1

% delay the sequence
max_num_delay_samps = ceil(D/Ts)-1;
y1 = delayseq([x5; zeros(max_num_delay_samps,1)], D, fs);
y1 = y1 + sig*(randn(length(y1),1)+1j*randn(length(y1),1));

% add fractional delayed multipaths
y5 = zeros(sps1, length(y1)+2);
for ii = 1:sps1
    y2 = resample(y1, sps1, 1);
    h = [1, zeros(1, ii-1), -0.95];
    hc{ii} = h;
    y3 = conv(y2,h);
    y4 = resample(y3,1,sps1);
    y5(ii,1:length(y4)) = y4;
end

% add integer delayed multipaths
int_max_delay = 11; % samples
step1 = 3;
y7 = zeros(ceil(int_max_delay/step1), length(y1)+int_max_delay);
kk = 1;
for ii = 1:step1:int_max_delay
    delays(kk) = ii;
    h = [1 zeros(1, ii-1) -0.95];
    hc2{kk} = h;
    y6 = conv(y1,h);
    y7(kk,1:length(y6)) = y6;
    kk = kk + 1;
end

%% Signal Shape Plots
figure
subplot(4,1,1)
plot(real(x1), 'b.-'); hold on
plot(imag(x1), 'rx-')
xlabel('Sample Number')
ylabel('Amplitude')
title('Input')
axis([0 inf -inf inf])

subplot(4,1,2)
plot(real(x5), 'b.-'); hold on
% plot(imag(x3), 'rx-')
xlabel('Sample Number')
ylabel('Amplitude')
title('Multipath Delayed Output')
axis([0 inf -inf inf])

subplot(4,1,3)
plot(real(y1), 'b.-')
xlabel('Sample Number')
ylabel('Amplitude')
title('Delayed Output')
axis([0 inf -inf inf])

subplot(4,1,4)
faxis = fs*(-0.5:1/1024:0.5-1/1024);
x5_f = 20*log10(fftshift(abs(fft(x5,1024))));
plot(faxis*1e-6, x5_f); hold on
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlabel('Frequency (MHz)')
ylabel('|fft|^2')
title('Received Signal Power Spectrum')


%% Plots for fractional delays
[z, lags] = xcorr(x5,x5);
[pval, pidx] = max(abs(z));
z_skew = skewness(abs(z(pidx-sps2:pidx+sps2)));

Nfft = max(1024, 2^(nextpow2(length(z)))); % must be even for correct results
faxis = fs*(-0.5:1/Nfft:0.5-1/Nfft);
bw = 0.9*fbw; %fs/sps2; % the bandwidth to use for computations
array_idx1 = find(faxis > -bw/2 & faxis < bw/2);
removal_band = 0.2e6; % band around zero to ignore, no delay information at k=0
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
H = [fk.' ones(length(fk),1)]; % used for LS estimation of delay

pdm = 0.15; % peak distance multipier
num_rows = sps1+2;
figure
subplot(num_rows,2,1)
plot(lags*Ts*1e9, abs(z), '.-'); hold on
plot(lags(pidx)*Ts*1e9, pval, 'rx')
plot(lags(pidx-sps2:pidx+sps2)*Ts*1e9, abs(z(pidx-sps2:pidx+sps2)), 'o', 'MarkerSize', 4)
axis([-pdm*N*sps2*Ts*1e9 pdm*N*sps2*Ts*1e9 0 1.2*max(abs(z))])
xlims = xlim;
xticks([xlims(1):Ts*1e9:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Fractional Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))
  
subplot(num_rows,2,2)
z_f = fft(x5,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), pi);

p = phase_f(array_idx);
m = H\p;
p_est = H*m;
delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
mean_delay_est = mean(delay_phase_est);

delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
plot(faxis*1e-6, phase_f, 'b-'); hold on
plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m(2), 'b.')
plot(faxis(array_idx)*1e-6, p_est, 'b-')
axis([faxis(1)*1e-6 faxis(end)*1e-6 -1 1])
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlims = xlim;
ylims = ylim;
text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f', mean_delay_est))
ylabel('Phase')
xlabel('Frequency (MHz)')

subplot(num_rows,2,3)
[z, lags] = xcorr(y1,x5);
[pval, pidx] = max(abs(z));
z_skew = skewness(abs(z(pidx-sps2:pidx+sps2)));
plot(lags*Ts*1e9, abs(z), '.-'); hold on
plot([D D]*1e9, ylim, 'k--')
plot(lags(pidx)*Ts*1e9, pval, 'rx')
plot(lags(pidx-sps2:pidx+sps2)*Ts*1e9, abs(z(pidx-sps2:pidx+sps2)), 'o', 'MarkerSize', 4)
axis([-pdm*N*sps2*Ts*1e9 pdm*N*sps2*Ts*1e9 0 1.2*max(abs(z))])
xlims = xlim;
xticks([xlims(1):Ts*1e9:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

subplot(num_rows,2,4)
% z_f = fft(ifftshift(z));
z_f = fft(y1,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), pi);

p = phase_f(array_idx);
m = H\p;
p_est = H*m-m(2);
delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
mean_delay_est = mean(delay_phase_est);

delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
plot(faxis*1e-6, phase_f-m(2), 'b-'); hold on
plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m(2), 'b.')
plot(faxis(array_idx)*1e-6, p_est, 'b-')
% axis([faxis(1)*1e-6 faxis(end)*1e-6 -5 5])
axis tight
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlims = xlim;ylims = ylim;
text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f', mean_delay_est))
ylabel('Phase')
xlabel('Frequency (MHz)')

kk = 1;
for ii = 1:2:2*(num_rows-2)
    subplot(num_rows,2,ii+4)
    [z, lags] = xcorr(y5(kk,:).',x5);
    [pval, pidx] = max(abs(z));
    z_skew = skewness(abs(z(pidx-sps2:pidx+sps2)));
    plot(lags*Ts*1e9, abs(z), '.-'); hold on
    plot([D D]*1e9, ylim, 'k--')
    plot(lags(pidx)*Ts*1e9, pval, 'rx')
    plot(lags(pidx-sps2:pidx+sps2)*Ts*1e9, abs(z(pidx-sps2:pidx+sps2)), 'o', 'MarkerSize', 4)
    axis([-pdm*N*sps2*Ts*1e9 pdm*N*sps2*Ts*1e9 0 1.2*max(abs(z))])
    xlims = xlim;
    xticks([xlims(1):Ts*1e9:xlims(2)])
    text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
    text(xlims(1), 0.75*pval, sprintf('Channel: h = [%s]', num2str(hc{kk},'%2.1f ')))
    xlabel('Lag (ns)')
    ylabel('|xcorr|')
    title(sprintf('Second Path Delay: %3.1f ns', kk*Ts/sps1*1e9))
    
%     z_f = fft(ifftshift(z));
    z_f = fft(y5(kk,:).',Nfft).*conj(fft(x5,Nfft));
    phase_f = unwrap(angle(fftshift(z_f)), pi);
    
    p = phase_f(array_idx);
    m = H\p;
    p_est = H*m-m(2);
    delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    mean_delay_est = mean(delay_phase_est);
    
    delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    subplot(num_rows,2,ii+5)
    plot(faxis*1e-6, phase_f-m(2), 'b-'); hold on
    plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m(2), 'b.')
    plot(faxis(array_idx)*1e-6, p_est, 'b-')
%     axis([faxis(1)*1e-6 faxis(end)*1e-6 -5 5])
    axis tight
    plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
    plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
    xlims = xlim;
    ylims = ylim;
    text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f', mean_delay_est))
    ylabel('Phase')
    xlabel('Frequency (MHz)')
        
    kk = kk + 1;
end

%% Plots for integer delays
[z, lags] = xcorr(x5,x5);
[pval, pidx] = max(abs(z));
z_skew = skewness(abs(z(pidx-sps2:pidx+sps2)));

Nfft = max(1024, 2^(nextpow2(length(z)))); % must be even for correct results
faxis = fs*(-0.5:1/Nfft:0.5-1/Nfft);
bw = 0.9*fbw; %fs/sps2; % the bandwidth to use for computations
array_idx1 = find(faxis > -bw/2 & faxis < bw/2);
removal_band = 0.2e6; % band around zero to ignore, no delay information at k=0
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
H = [fk.' ones(length(fk),1)]; % used for LS estimation of delay

pdm = 0.15; % peak distance multipier
num_rows = ceil(int_max_delay/step1)+2;
figure
subplot(num_rows,2,1)
plot(lags*Ts*1e9, abs(z), '.-'); hold on
plot(lags(pidx)*Ts*1e9, pval, 'rx')
plot(lags(pidx-sps2:pidx+sps2)*Ts*1e9, abs(z(pidx-sps2:pidx+sps2)), 'o', 'MarkerSize', 4)
axis([-pdm*N*sps2*Ts*1e9 pdm*N*sps2*Ts*1e9 0 1.2*max(abs(z))])
xlims = xlim;
xticks([xlims(1):Ts*1e9:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Integer Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))
  
subplot(num_rows,2,2)
z_f = fft(x5,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), pi);

p = phase_f(array_idx);
m = H\p;
p_est = H*m;
delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
mean_delay_est = mean(delay_phase_est);

delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
plot(faxis*1e-6, phase_f, 'b-'); hold on
plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m(2), 'b.')
plot(faxis(array_idx)*1e-6, p_est, 'b-')
axis([faxis(1)*1e-6 faxis(end)*1e-6 -1 1])
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlims = xlim;
ylims = ylim;
text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f', mean_delay_est))
ylabel('Phase')
xlabel('Frequency (MHz)')

subplot(num_rows,2,3)
[z, lags] = xcorr(y1,x5);
[pval, pidx] = max(abs(z));
z_skew = skewness(abs(z(pidx-sps2:pidx+sps2)));
plot(lags*Ts*1e9, abs(z), '.-'); hold on
plot([D D]*1e9, ylim, 'k--')
plot(lags(pidx)*Ts*1e9, pval, 'rx')
plot(lags(pidx-sps2:pidx+sps2)*Ts*1e9, abs(z(pidx-sps2:pidx+sps2)), 'o', 'MarkerSize', 4)
axis([-pdm*N*sps2*Ts*1e9 pdm*N*sps2*Ts*1e9 0 1.2*max(abs(z))])
xlims = xlim;
xticks([xlims(1):Ts*1e9:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

subplot(num_rows,2,4)
% z_f = fft(ifftshift(z));
z_f = fft(y1,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), pi);

p = phase_f(array_idx);
m = H\p;
p_est = H*m-m(2);
delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
mean_delay_est = mean(delay_phase_est);

delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
plot(faxis*1e-6, phase_f-m(2), 'b-'); hold on
plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m(2), 'b.')
plot(faxis(array_idx)*1e-6, p_est, 'b-')
% axis([faxis(1)*1e-6 faxis(end)*1e-6 -5 5])
axis tight
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlims = xlim;ylims = ylim;
text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f', mean_delay_est))
ylabel('Phase')
xlabel('Frequency (MHz)')

kk = 1;
for ii = 1:2:2*(num_rows-2)
    subplot(num_rows,2,ii+4)
    [z, lags] = xcorr(y7(kk,:).',x5);
    [pval, pidx] = max(abs(z));
    z_skew = skewness(abs(z(pidx-sps2:pidx+sps2)));
    plot(lags*Ts*1e9, abs(z), '.-'); hold on
    plot([D D]*1e9, ylim, 'k--')
    plot(lags(pidx)*Ts*1e9, pval, 'rx')
    plot(lags(pidx-sps2:pidx+sps2)*Ts*1e9, abs(z(pidx-sps2:pidx+sps2)), 'o', 'MarkerSize', 4)
    axis([-pdm*N*sps2*Ts*1e9 pdm*N*sps2*Ts*1e9 0 1.2*max(abs(z))])
    xlims = xlim;
    xticks([xlims(1):Ts*1e9:xlims(2)])
    text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
    text(xlims(1), 0.75*pval, sprintf('Channel: h = [%s]', num2str(hc{kk},'%2.1f ')))
    xlabel('Lag (ns)')
    ylabel('|xcorr|')
    title(sprintf('Second Path Delay: %3.1f ns', delays(kk)))
    
%     z_f = fft(ifftshift(z));
    z_f = fft(y7(kk,:).',Nfft).*conj(fft(x5,Nfft));
    phase_f = unwrap(angle(fftshift(z_f)), pi);
    
    p = phase_f(array_idx);
    m = H\p;
    p_est = H*m-m(2);
    delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    mean_delay_est = mean(delay_phase_est);
    
    delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    subplot(num_rows,2,ii+5)
    plot(faxis*1e-6, phase_f-m(2), 'b-'); hold on
    plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m(2), 'b.')
    plot(faxis(array_idx)*1e-6, p_est, 'b-')
%     axis([faxis(1)*1e-6 faxis(end)*1e-6 -5 5])
    axis tight
    plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
    plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
    xlims = xlim;
    ylims = ylim;
    text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f', mean_delay_est))
    ylabel('Phase')
    xlabel('Frequency (MHz)')
        
    kk = kk + 1;
end



% %% Plots for integer delays - old
% figure
% pdm = 0.35; % peak distance multipier
% num_rows = length(delays)+2;
% subplot(num_rows,2,1)
% [z, lags] = xcorr(x5,x5);
% [pval, pidx] = max(abs(z));
% z_skew = skewness(abs(z(pidx-sps2:pidx+sps2)));
% plot(lags*Ts*1e9, abs(z), '.-'); hold on
% plot(lags(pidx)*Ts*1e9, pval, 'rx')
% plot(lags(pidx-sps2:pidx+sps2)*Ts*1e9, abs(z(pidx-sps2:pidx+sps2)), 'o', 'MarkerSize', 4)
% axis([-pdm*N*sps2*Ts*1e9 pdm*N*sps2*Ts*1e9 0 1.2*max(abs(z))])
% xlims = xlim;
% xticks([xlims(1):Ts*1e9:xlims(2)])
% text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
% xlabel('Lag (ns)')
% ylabel('|xcorr|')
% title(sprintf('Integer Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))
% 
% subplot(num_rows,2,2)
% z_f = fft(z, Nfft);
% phase_f = unwrap(angle(fftshift(z_f)), pi);
% delay_phase = -Nfft*phase_f(:)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
% plot(faxis*1e-6, phase_f); hold on
% plot(faxis*1e-6, delay_phase, '.')
% plot([-fs/2 -fs/2]/sps2*1e-6, ylim, 'k--')
% plot([fs/2 fs/2]/sps2*1e-6, ylim, 'k--')
% ylabel('Phase')
% xlabel('Frequency (MHz)')
% axis([-inf inf -1000 300])
%     
% subplot(num_rows,2,3)
% [z, lags] = xcorr(y1,x5);
% [pval, pidx] = max(abs(z));
% z_skew = skewness(abs(z(pidx-sps2:pidx+sps2)));
% plot(lags*Ts*1e9, abs(z), '.-'); hold on
% plot([D D]*1e9, ylim, 'k--')
% plot(lags(pidx)*Ts*1e9, pval, 'rx')
% plot(lags(pidx-sps2:pidx+sps2)*Ts*1e9, abs(z(pidx-sps2:pidx+sps2)), 'o', 'MarkerSize', 4)
% axis([-pdm*N*sps2*Ts*1e9 pdm*N*sps2*Ts*1e9 0 1.2*max(abs(z))])
% xlims = xlim;
% xticks([xlims(1):Ts*1e9:xlims(2)])
% text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
% xlabel('Lag (ns)')
% ylabel('|xcorr|')
% title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))
% 
% subplot(num_rows,2,4)
% z_f = fft(z, Nfft);
% phase_f = unwrap(angle(fftshift(z_f)));
% delay_phase = -Nfft*phase_f(:)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
% plot(faxis*1e-6, phase_f); hold on
% plot(faxis*1e-6, delay_phase, '.')
% plot([-fs/2 -fs/2]/sps2*1e-6, ylim, 'k--')
% plot([fs/2 fs/2]/sps2*1e-6, ylim, 'k--')
% ylabel('Phase')
% xlabel('Frequency (MHz)')
% axis([-inf inf -1000 300])
% 
% kk = 1;
% for ii = 1:2:2*(num_rows-2)
%     subplot(num_rows,2,ii+4)
%     [z, lags] = xcorr(y7(kk,:),x5);
%     [pval, pidx] = max(abs(z));
%     z_skew = skewness(abs(z(pidx-sps2:pidx+sps2)),0);
%     plot(lags*Ts*1e9, abs(z), '.-'); hold on
%     plot([D D]*1e9, ylim, 'k--')
%     plot(lags(pidx)*Ts*1e9, pval, 'rx')
%     plot(lags(pidx-sps2:pidx+sps2)*Ts*1e9, abs(z(pidx-sps2:pidx+sps2)), 'o', 'MarkerSize', 4)
%     axis([-pdm*N*sps2*Ts*1e9 pdm*N*sps2*Ts*1e9 0 1.2*max(abs(z))])
%     xlims = xlim;
%     xticks([xlims(1):200:-1200 -1000:Ts*1e9:1000 1200:200:xlims(2)])
%     text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
%     text(xlims(1), 0.75*pval, sprintf('Channel: h = [%s]', num2str(hc2{kk},'%2.1f ')))
%     xlabel('Lag (ns)')
%     ylabel('|xcorr|')
%     title(sprintf('Second Path Delay: %3.1f ns', delays(kk)*Ts*1e9))
%     
%     z_f = fft(z, Nfft);
%     phase_f = unwrap(angle(fftshift(z_f)));
%     delay_phase = -Nfft*phase_f(:)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
%     subplot(num_rows,2,ii+5)
%     plot(faxis*1e-6, phase_f); hold on
%     plot(faxis*1e-6, delay_phase, '.')
%     plot([-fs/2 -fs/2]/sps2*1e-6, ylim, 'k--')
%     plot([fs/2 fs/2]/sps2*1e-6, ylim, 'k--')
%     ylabel('Phase')
%     xlabel('Frequency (MHz)')
%     axis([-inf inf -1000 300])
%     
%     kk = kk + 1;
% end


function [taps] = exponential(Ntaps, delay_spread)

f = @(x) exp(-x/delay_spread);

taps = f(0:Ntaps-1);
phi = pi*(2*rand(1,Ntaps)-1);
taps = taps.*exp(1j*phi);

end
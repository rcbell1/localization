clear; close all

N = 100;     % length of first sequence (symbols)
sps = 4;   % samples per symbol for shaping filter
fbw = 7e6;  % occupied bandwidth of transmitter
% fs = fbw*sps; % 14e6;  % receiver sampling rate (Hz)
fs = 28e6;  % receiver sampling rate (Hz)
Ts = 1/fs;  % sample period (s)
D = 3*(2*rand-1)*Ts;  % delay (s)
sig = 0.0005; % standard deviation of noise
multi_sps = 3;  % oversample factor when multipath is added
% sps = 4;   % samples per symbol for shaping filter
M = 2;      % modulation order
pulse_shaping = 0;  % pulse shape using root raised cosine 0/1 - no/yes
span = 5;   % length of pulse shaping filter (symbols)
beta = 0.18; % excess bandwidth of pulse shaping filter

% Cyclostationary plot parameters (normalized units)
os_fac    = 4; % oversampling amount
os_fac_conj = 4;
sym_rate = fbw/fs; % noramalized symbol rate
alpha_scd = 0.:sym_rate/4:1; % define cycle frequencies to compute
alpa_sym_idx = find(alpha_scd == sym_rate);
% alpha_scd = [0 0.3 0.5 0.7];
freq_res = 0.15;
wind_o_n  = os_fac*round(freq_res*N*sps); % number of window taps
wind_o    = ones(wind_o_n,1)/wind_o_n; % normalized window

% create the baseline modulated signal
xb = randi([0 M-1], log2(M)*N, 1);
x1 = qammod(xb, M, 'gray');
x2 = upsample(x1, sps);
if pulse_shaping == 1
    rrc = rcosdesign(beta, span, sps, 'sqrt');
    shape_filt = rrc.'/max(rrc);
else
    shape_filt = ones(sps,1);
end
x3 = filter(shape_filt, 1, x2);

[P,Q] = rat((fs/fbw)/sps);
if P*Q == 1
    tx_sig = x3/sqrt(mean(abs(x3).^2)); % normalize power to 1
else
    x4 = resample(x3, P, Q);
    tx_sig = x4/sqrt(mean(abs(x4).^2)); % normalize power to 1
end

% x4=x3; % debug
% x4 = sin(2*pi*1/sps2*(0:sps2*N-1)).'; % debug

% delay the sequence
max_num_delay_samps = ceil(D/Ts)-1;
rx_sig_no_multi = delayseq([tx_sig; zeros(max_num_delay_samps,1)], D, fs);
rx_sig_no_multi = rx_sig_no_multi + sig*(randn(length(rx_sig_no_multi),1) + ...
    1j*randn(length(rx_sig_no_multi),1));

% add fractional delayed multipaths
multi_amp = 0.3;
[rx_sig_frac_multi, int_channels] = add_frac_multi(rx_sig_no_multi, multi_sps, multi_amp);

% add integer delayed multipaths
int_max_delay = 11; % samples
step1 = 4;
[rx_sig_int_multi, delays, frac_channels] = ...
    add_integer_multi(rx_sig_no_multi, int_max_delay, step1, multi_amp);

%% Signal Shape Plots
figure
subplot(4,1,1)
plot(real(x1), 'b.-'); hold on
plot(imag(x1), 'r.-')
xlabel('Sample Number')
ylabel('Amplitude')
title('Input')
axis([0 inf -inf inf])

subplot(4,1,2)
% plot(real(x5), 'b.-'); hold on
plot(real(x3), 'b.-')
xlabel('Sample Number')
ylabel('Amplitude')
title('Shaped Channel Input')
axis([0 inf -inf inf])

subplot(4,1,3)
plot(real(rx_sig_no_multi), 'b.-')
xlabel('Sample Number')
ylabel('Amplitude')
title('Delayed AWGN Channel Output')
axis([0 inf -inf inf])

subplot(4,1,4)
Nfft = max(1024, 2^(nextpow2(size(rx_sig_int_multi,2)))); % must be even for correct results
faxis = fs*(-0.5:1/Nfft:0.5-1/Nfft);
tx_sig_f = 10*log10(fftshift(1/Nfft*abs(fft(tx_sig,Nfft)).^2));
rx_sig_no_multi_f = 20*log10(fftshift(1/Nfft*abs(fft(rx_sig_no_multi.',Nfft)).^2));
rx_sig_frac_multi_f = 20*log10(fftshift(1/Nfft*abs(fft(rx_sig_frac_multi.',Nfft)).^2));
rx_sig_int_multi_f = 20*log10(fftshift(1/Nfft*abs(fft(rx_sig_int_multi.',Nfft)).^2));
lh = [];
lh(1) = plot(faxis*1e-6, tx_sig_f, '.-'); hold on
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlabel('Frequency (MHz)')
ylabel('|fft|^2')
title('Signal Power Spectrum')

%% Cyclostationarity check of input sequence
nfft      = 2^(nextpow2(os_fac*length(x3)));
norm_f_vec = (1:nfft) - nfft/2;
norm_f_vec = norm_f_vec/max(2*abs(norm_f_vec));

scd_fsm_mat = os_fac*cyclo_scd_fsm(tx_sig,alpha_scd,wind_o,nfft);
% scd_conj_fsm_mat = os_fac_conj*cyclo_scd_fsm_conj(tx_sig,alpha_scd,wind_o,nfft);

% 3D plot of spectral freq by cycle freq by SCD magnitude
Z = abs(scd_fsm_mat);
plot3d_slice(norm_f_vec, alpha_scd, Z, sym_rate);
text(0.5,sym_rate,0.5,sprintf('Symbol Rate: %.2f',sym_rate), 'Rotation',-25, 'VerticalAlignment', 'bottom')

Z = angle(scd_fsm_mat);
Z = Z([1,alpa_sym_idx],:); % alpha = 0 and sym_rate
plot3d_slice(norm_f_vec, alpha_scd([1,alpa_sym_idx]), Z, sym_rate);
text(0.5,0,0.5,sprintf('Symbol Rate: %.2f',sym_rate), 'Rotation',-25, 'VerticalAlignment', 'bottom')


%% Plots for fractional delays
[z, lags] = xcorr(tx_sig,tx_sig);
[pval, pidx] = max(abs(z));
z_skew = skewness(abs(z(pidx-sps:pidx+sps)));

Nfft = max(1024, 2^(nextpow2(length(z)))); % must be even for correct results
faxis = fs*(-0.5:1/Nfft:0.5-1/Nfft);
tx_sig_f = 10*log10(fftshift(1/Nfft*abs(fft(tx_sig,Nfft)).^2));
rx_sig_no_multi_f = 20*log10(fftshift(1/Nfft*abs(fft(rx_sig_no_multi.',Nfft)).^2));
rx_sig_frac_multi_f = 20*log10(fftshift(1/Nfft*abs(fft(rx_sig_frac_multi.',Nfft)).^2));
rx_sig_int_multi_f = 20*log10(fftshift(1/Nfft*abs(fft(rx_sig_int_multi.',Nfft)).^2));

% estimate the delay using phase
bw = 0.9*fbw; %fs/sps2; % the bandwidth to use for computations
removal_band = 0.2e6; % band around zero to ignore, no delay information at k=0
phase_wrap_step = pi;

[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(tx_sig, tx_sig, faxis, Nfft, fs, bw, removal_band, phase_wrap_step);

% array_idx2 = find(faxis > -removal_band & faxis < removal_band);
% [~, shared_idx] = intersect(array_idx1,array_idx2,'stable');
% array_idx1(shared_idx) = []; % delete those frequency indices in the removal band
% array_idx = array_idx1;
% kidx = faxis*Nfft/fs;
% fkeep_idx = find(faxis > -bw/2 & faxis < bw/2);
% frem_idx = find(faxis > -removal_band & faxis < removal_band);
% [~, shared_f] = intersect(fkeep_idx,frem_idx,'stable');
% fkeep_idx(shared_idx) = []; % delete those frequency indices in the removal band
% fk = faxis(fkeep_idx);
% kidx = kidx(fkeep_idx);
% 
% z_f = fft(tx_sig,Nfft).*conj(fft(tx_sig,Nfft));
% phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);
% H = [fk.' ones(length(fk),1)]; % used for LS estimation of delay
% p = phase_f(array_idx);
% m = H\p;
% p_est = H*m;
% delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
% mean_delay_est = mean(delay_phase_est);



pdm = 0.15; % peak distance multipier just for plot display
xtick_step = 4*Ts*1e9;
num_rows = multi_sps+2;
figure
subplot(num_rows,2,1)
plot(lags*Ts*1e9, abs(z), '.-'); hold on
plot(lags(pidx)*Ts*1e9, pval, 'rx')
plot(lags(pidx-sps:pidx+sps)*Ts*1e9, abs(z(pidx-sps:pidx+sps)), 'o', 'MarkerSize', 4)
axis([-pdm*N*sps*Ts*1e9 pdm*N*sps*Ts*1e9 0 1.2*max(abs(z))])
xlims = xlim;
xticks([xlims(1):xtick_step:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Fractional Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))
  
subplot(num_rows,2,2)
% z_f = fft(tx_sig,Nfft).*conj(fft(tx_sig,Nfft));
% phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);

% H = [fk.' ones(length(fk),1)]; % used for LS estimation of delay
% p = phase_f(array_idx);
% m = H\p;
% p_est = H*m;
% delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
% mean_delay_est = mean(delay_phase_est);

% delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
plot(faxis*1e-6, phase_f, 'b-'); hold on
plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m, 'b.')
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
[z, lags] = xcorr(rx_sig_no_multi,tx_sig);
[pval, pidx] = max(abs(z));
z_skew = skewness(abs(z(pidx-sps:pidx+sps)));

plot(lags*Ts*1e9, abs(z), '.-'); hold on
plot([D D]*1e9, ylim, 'k--')
plot(lags(pidx)*Ts*1e9, pval, 'rx')
plot(lags(pidx-sps:pidx+sps)*Ts*1e9, abs(z(pidx-sps:pidx+sps)), 'o', 'MarkerSize', 4)
axis([-pdm*N*sps*Ts*1e9 pdm*N*sps*Ts*1e9 0 1.2*max(abs(z))])
xlims = xlim;
xticks([xlims(1):xtick_step:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

subplot(num_rows,2,4)
% z_f = fft(ifftshift(z));
% z_f = fft(rx_sig_no_multi,Nfft).*conj(fft(tx_sig,Nfft));
% phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);

[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_no_multi, tx_sig, faxis, Nfft, fs, bw, removal_band, phase_wrap_step);

% p = phase_f(array_idx);
% m = H\p;
% p_est = H*m-m(2);
% delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
% mean_delay_est = mean(delay_phase_est);

% delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
plot(faxis*1e-6, phase_f-m, 'b-'); hold on
plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m, 'b.')
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
    % plot correlation in time
    subplot(num_rows,2,ii+4)
    [z, lags] = xcorr(rx_sig_frac_multi(kk,:).',tx_sig);
    [pval, pidx] = max(abs(z));
    z_skew = skewness(abs(z(pidx-sps:pidx+sps)));
    plot(lags*Ts*1e9, abs(z), '.-'); hold on
    plot([D D]*1e9, ylim, 'k--')
    plot(lags(pidx)*Ts*1e9, pval, 'rx')
    plot(lags(pidx-sps:pidx+sps)*Ts*1e9, abs(z(pidx-sps:pidx+sps)), 'o', 'MarkerSize', 4)
    axis([-pdm*N*sps*Ts*1e9 pdm*N*sps*Ts*1e9 0 1.2*max(abs(z))])
    xlims = xlim;
    xticks([xlims(1):xtick_step:xlims(2)])
    text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
    text(xlims(1), 0.75*pval, sprintf('Channel: h = [%s]', num2str(frac_channels{kk},'%2.1f ')))
    xlabel('Lag (ns)')
    ylabel('|xcorr|')
    title(sprintf('Second Path Delay: %3.1f ns', kk*Ts/multi_sps*1e9))
    
    % plot phase of correlation in frequency
    subplot(num_rows,2,ii+5)
    [mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_frac_multi(kk,:).', tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);
%     z_f = fft(ifftshift(z));
%     z_f = fft(rx_sig_frac_multi(kk,:).',Nfft).*conj(fft(tx_sig,Nfft));
%     phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);
%     
%     p = phase_f(array_idx);
%     m = H\p;
%     p_est = H*m-m(2);
%     delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
%     mean_delay_est = mean(delay_phase_est);
%     
%     delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)

    plot(faxis*1e-6, phase_f-m, 'b-'); hold on
    plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m, 'b.')
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

%% Fractional Delay Magnitude ontop of Phase Plots
pdm = 0.15; % peak distance multipier for plots
% phase_wrap_step = pi;
xtick_step = 4*Ts*1e9;
num_rows = multi_sps+2;

figure
subplot(num_rows,1,1)
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(tx_sig, tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);  
% z_f = fft(tx_sig,Nfft).*conj(fft(tx_sig,Nfft));
% phase_f = unwrap(angle(fftshift(z_f)), pi);
% 
% p = phase_f(array_idx);
% m = H\p;
% p_est = H*m;
% delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
% mean_delay_est = mean(delay_phase_est);
% 
% delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
yyaxis left
plot(faxis*1e-6, tx_sig_f, 'b-'); hold on
axis([-inf inf -100 0])
ylabel('PSD')
yyaxis right
plot(faxis*1e-6, phase_f, 'r-')
axis([faxis(1)*1e-6 faxis(end)*1e-6 -1 1])
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlims = xlim;
ylims = ylim;
ylabel('Phase')
xlabel('Frequency (MHz)')
title(sprintf('Fractional Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))

subplot(num_rows,1,2)
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_no_multi, tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);    
% z_f = fft(rx_sig_no_multi,Nfft).*conj(fft(tx_sig,Nfft));
% phase_f = unwrap(angle(fftshift(z_f)), pi);
% 
% p = phase_f(array_idx);
% m = H\p;
% p_est = H*m-m(2);
% delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
% mean_delay_est = mean(delay_phase_est);
% 
% delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
yyaxis left
plot(faxis*1e-6, rx_sig_no_multi_f, 'b-'); hold on
axis([-inf inf -100 0])
ylabel('PSD')
yyaxis right
plot(faxis*1e-6, phase_f-m, 'r-')
axis tight
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlims = xlim;ylims = ylim;
ylabel('Phase')
xlabel('Frequency (MHz)')
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

kk = 1;
for ii = 1:1:1*(num_rows-2)
    subplot(num_rows,1,ii+2)   
    [mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_frac_multi(kk,:).', tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);   
%     z_f = fft(rx_sig_frac_multi(kk,:).',Nfft).*conj(fft(tx_sig,Nfft));
%     phase_f = unwrap(angle(fftshift(z_f)), pi);
%     
%     p = phase_f(array_idx);
%     m = H\p;
%     p_est = H*m-m(2);
%     delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
%     mean_delay_est = mean(delay_phase_est);
%     
%     delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    yyaxis left
    plot(faxis*1e-6, rx_sig_frac_multi_f(:,ii), 'b-'); hold on
    axis([-inf inf -100 0])
    ylabel('PSD')
    yyaxis right
    plot(faxis*1e-6, phase_f-m, 'r-')
    axis tight
    plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
    plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
    xlims = xlim;
    ylims = ylim;
    ylabel('Phase')
    xlabel('Frequency (MHz)')
%     title(sprintf('Delay Spread: %3.1f ns', delspreads(kk)*1e9))
        
    kk = kk + 1;
end

%% ***************** Cyclostationarity plots *************************
% *******************************************************************
scd_fsm_mat = os_fac*cyclo_cross_scd_fsm(tx_sig,rx_sig_no_multi,alpha_scd,wind_o,nfft);
Z = abs(scd_fsm_mat);
plot3d_slice(norm_f_vec, alpha_scd, Z, sym_rate);
text(0.5,sym_rate,0.5,sprintf('Symbol Rate: %.2f',sym_rate), 'Rotation',-25, 'VerticalAlignment', 'bottom')

Z = angle(scd_fsm_mat);
Z = Z([1,alpa_sym_idx],:); % alpha = 0 and sym_rate
plot3d_slice(norm_f_vec, alpha_scd([1,5]), Z, sym_rate);
text(0.5,0,0.5,sprintf('Symbol Rate: %.2f',sym_rate), 'Rotation',-25, 'VerticalAlignment', 'bottom')

%% Plots for integer delays
[z, lags] = xcorr(tx_sig,tx_sig);
[pval, pidx] = max(abs(z));
z_skew = skewness(abs(z(pidx-sps:pidx+sps)));

Nfft = max(1024, 2^(nextpow2(length(z)))); % must be even for correct results
faxis = fs*(-0.5:1/Nfft:0.5-1/Nfft);
% bw = 0.9*fbw; %fs/sps2; % the bandwidth to use for computations
% array_idx1 = find(faxis > -bw/2 & faxis < bw/2);
% removal_band = 0.2e6; % band around zero to ignore, no delay information at k=0
% array_idx2 = find(faxis > -removal_band & faxis < removal_band);
% [~, shared_idx] = intersect(array_idx1,array_idx2,'stable');
% array_idx1(shared_idx) = []; % delete those frequency indices in the removal band
% array_idx = array_idx1;
% kidx = faxis*Nfft/fs;
% fkeep_idx = find(faxis > -bw/2 & faxis < bw/2);
% frem_idx = find(faxis > -removal_band & faxis < removal_band);
% [~, shared_f] = intersect(fkeep_idx,frem_idx,'stable');
% fkeep_idx(shared_idx) = []; % delete those frequency indices in the removal band
% fk = faxis(fkeep_idx);
% kidx = kidx(fkeep_idx);
% H = [fk.' ones(length(fk),1)]; % used for LS estimation of delay

pdm = 0.15; % peak distance multipier
num_rows = ceil(int_max_delay/step1)+2;
figure
subplot(num_rows,2,1)
plot(lags*Ts*1e9, abs(z), '.-'); hold on
plot(lags(pidx)*Ts*1e9, pval, 'rx')
plot(lags(pidx-sps:pidx+sps)*Ts*1e9, abs(z(pidx-sps:pidx+sps)), 'o', 'MarkerSize', 4)
axis([-pdm*N*sps*Ts*1e9 pdm*N*sps*Ts*1e9 0 1.2*max(abs(z))])
xlims = xlim;
xticks([xlims(1):xtick_step:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Integer Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))
  
subplot(num_rows,2,2)
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(tx_sig, tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);  
% z_f = fft(tx_sig,Nfft).*conj(fft(tx_sig,Nfft));
% phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);
% 
% p = phase_f(array_idx);
% m = H\p;
% p_est = H*m;
% delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
% mean_delay_est = mean(delay_phase_est);
% 
% delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
plot(faxis*1e-6, phase_f, 'b-'); hold on
plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m, 'b.')
plot(faxis(array_idx)*1e-6, p_est, 'b-')
axis([faxis(1)*1e-6 faxis(end)*1e-6 -1 1])
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlims = xlim;
ylims = ylim;
text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f ns', mean_delay_est))
ylabel('Phase')
xlabel('Frequency (MHz)')

subplot(num_rows,2,3)
[z, lags] = xcorr(rx_sig_no_multi,tx_sig);
[pval, pidx] = max(abs(z));
z_skew = skewness(abs(z(pidx-sps:pidx+sps)));
plot(lags*Ts*1e9, abs(z), '.-'); hold on
plot([D D]*1e9, ylim, 'k--')
plot(lags(pidx)*Ts*1e9, pval, 'rx')
plot(lags(pidx-sps:pidx+sps)*Ts*1e9, abs(z(pidx-sps:pidx+sps)), 'o', 'MarkerSize', 4)
axis([-pdm*N*sps*Ts*1e9 pdm*N*sps*Ts*1e9 0 1.2*max(abs(z))])
xlims = xlim;
xticks([xlims(1):xtick_step:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

subplot(num_rows,2,4)
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_no_multi, tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);  
% z_f = fft(ifftshift(z));
% z_f = fft(rx_sig_no_multi,Nfft).*conj(fft(tx_sig,Nfft));
% phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);
% 
% p = phase_f(array_idx);
% m = H\p;
% p_est = H*m-m(2);
% delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
% mean_delay_est = mean(delay_phase_est);
% 
% delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
plot(faxis*1e-6, phase_f-m, 'b-'); hold on
plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m, 'b.')
plot(faxis(array_idx)*1e-6, p_est, 'b-')
% axis([faxis(1)*1e-6 faxis(end)*1e-6 -5 5])
axis tight
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlims = xlim;ylims = ylim;
text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f ns', mean_delay_est))
ylabel('Phase')
xlabel('Frequency (MHz)')

kk = 1;
for ii = 1:2:2*(num_rows-2)
    subplot(num_rows,2,ii+4)
    [z, lags] = xcorr(rx_sig_int_multi(kk,:).',tx_sig);
    [pval, pidx] = max(abs(z));
    z_skew = skewness(abs(z(pidx-sps:pidx+sps)));
    plot(lags*Ts*1e9, abs(z), '.-'); hold on
    plot([D D]*1e9, ylim, 'k--')
    plot(lags(pidx)*Ts*1e9, pval, 'rx')
    plot(lags(pidx-sps:pidx+sps)*Ts*1e9, abs(z(pidx-sps:pidx+sps)), 'o', 'MarkerSize', 4)
    axis([-pdm*N*sps*Ts*1e9 pdm*N*sps*Ts*1e9 0 1.2*max(abs(z))])
    xlims = xlim;
    xticks([xlims(1):xtick_step:xlims(2)])
    text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
    text(xlims(1), 0.75*pval, sprintf('Channel: h = [%s]', num2str(int_channels{kk},'%2.1f ')))
    xlabel('Lag (ns)')
    ylabel('|xcorr|')
    title(sprintf('Second Path Delay: %3.1f ns', delays(kk)*Ts*1e9))
    
%     z_f = fft(ifftshift(z));
%     z_f = fft(rx_sig_int_multi(kk,:).',Nfft).*conj(fft(tx_sig,Nfft));
%     phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);
%     
%     p = phase_f(array_idx);
%     m = H\p;
%     p_est = H*m-m(2);
%     delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
%     mean_delay_est = mean(delay_phase_est);
%     
%     delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    subplot(num_rows,2,ii+5)
    [mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_int_multi(kk,:).', tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);  

    plot(faxis*1e-6, phase_f-m, 'b-'); hold on
    plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m, 'b.')
    plot(faxis(array_idx)*1e-6, p_est, 'b-')
%     axis([faxis(1)*1e-6 faxis(end)*1e-6 -5 5])
    axis tight
    plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
    plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
    xlims = xlim;
    ylims = ylim;
    text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f ns', mean_delay_est))
    ylabel('Phase')
    xlabel('Frequency (MHz)')
        
    kk = kk + 1;
end

%% Integer Delay Magnitude ontop of Phase Plots
pdm = 0.15; % peak distance multipier for plots
% phase_wrap_step = pi;
xtick_step = 4*Ts*1e9;
num_rows = multi_sps+2;

figure
subplot(num_rows,1,1)
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(tx_sig, tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);  
% z_f = fft(tx_sig,Nfft).*conj(fft(tx_sig,Nfft));
% phase_f = unwrap(angle(fftshift(z_f)), pi);
% 
% p = phase_f(array_idx);
% m = H\p;
% p_est = H*m;
% delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
% mean_delay_est = mean(delay_phase_est);
% 
% delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
yyaxis left
plot(faxis*1e-6, tx_sig_f(:,1), 'b-'); hold on
ylabel('PSD')
axis([-inf inf -100 0])
yyaxis right
plot(faxis*1e-6, phase_f, 'r-')
axis([faxis(1)*1e-6 faxis(end)*1e-6 -1 1])
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlims = xlim;
ylims = ylim;
ylabel('Phase')
xlabel('Frequency (MHz)')
title(sprintf('Integer Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))

subplot(num_rows,1,2)
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_no_multi, tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);  
% z_f = fft(rx_sig_no_multi,Nfft).*conj(fft(tx_sig,Nfft));
% phase_f = unwrap(angle(fftshift(z_f)), pi);
% 
% p = phase_f(array_idx);
% m = H\p;
% p_est = H*m-m(2);
% delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
% mean_delay_est = mean(delay_phase_est);
% 
% delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
yyaxis left
plot(faxis*1e-6, rx_sig_no_multi_f, 'b-'); hold on
axis([-inf inf -100 0])
ylabel('PSD')
yyaxis right
plot(faxis*1e-6, phase_f-m, 'r-')
axis tight
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlims = xlim;ylims = ylim;
ylabel('Phase')
xlabel('Frequency (MHz)')
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

kk = 1;
for ii = 1:1:1*(num_rows-2)
    subplot(num_rows,1,ii+2)    
    [mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
        get_phase_est(rx_sig_int_multi(kk,:).', tx_sig, faxis, Nfft, fs, ...
        bw, removal_band, phase_wrap_step); 
%     z_f = fft(rx_sig_int_multi(kk,:).',Nfft).*conj(fft(tx_sig,Nfft));
%     phase_f = unwrap(angle(fftshift(z_f)), pi);
%     
%     p = phase_f(array_idx);
%     m = H\p;
%     p_est = H*m-m(2);
%     delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
%     mean_delay_est = mean(delay_phase_est);
%     
%     delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    yyaxis left
    plot(faxis*1e-6, rx_sig_int_multi_f(:,ii), 'b-'); hold on
    axis([-inf inf -100 0])
    ylabel('PSD')
    yyaxis right
    plot(faxis*1e-6, phase_f-m, 'r-')
    axis tight
    plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
    plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
    xlims = xlim;
    ylims = ylim;
    ylabel('Phase')
    xlabel('Frequency (MHz)')
    
    kk = kk + 1;
end

function [y_frac_multi, frac_channels] = add_frac_multi(in, frac_sps, multi_amp)
    y_frac_multi = zeros(frac_sps, length(in)+2);
    for ii = 1:frac_sps
        y2 = resample(in, frac_sps, 1);
        h = [1, zeros(1, ii-1), multi_amp];
        frac_channels{ii} = h;
        y3 = conv(y2,h);
        y4 = resample(y3,1,frac_sps);
        y_frac_multi(ii,1:length(y4)) = y4;
    end
end

function [y_int_multi, delays, int_channels] = add_integer_multi(in, int_max_delay, step, multi_amp)
    y_int_multi = zeros(ceil(int_max_delay/step), length(in)+int_max_delay);
    kk = 1;
    for ii = 1:step:int_max_delay
        delays(kk) = ii;
        h = [1 zeros(1, ii-1) multi_amp];
        int_channels{kk} = h;
        y2 = conv(in,h);
        y_int_multi(kk,1:length(y2)) = y2;
        kk = kk + 1;
    end
end

function [mean_delay_est, phase_f, ls_phase_est, ls_mean_phase_est, ...
    fkeep_freqs, occbw_idx, kidx] = ...
    get_phase_est(sig1, sig2, faxis, Nfft, fs, bw, removal_band, phase_wrap_step)

    occbw_idx = find(faxis > -bw/2 & faxis < bw/2);    
    removal_idx = find(faxis > -removal_band & faxis < removal_band);
    [~, shared_idx] = intersect(occbw_idx,removal_idx,'stable');
    occbw_idx(shared_idx) = []; % delete those frequency indices in the removal band
%     array_idx = occbw_idx;
    kidx = faxis*Nfft/fs;
    fkeep_idx = find(faxis > -bw/2 & faxis < bw/2);
    frem_idx = find(faxis > -removal_band & faxis < removal_band);
%     [~, shared_f] = intersect(fkeep_idx,frem_idx,'stable');
    fkeep_idx(shared_idx) = []; % delete those frequency indices in the removal band
    fkeep_freqs = faxis(fkeep_idx);
    kidx = kidx(fkeep_idx);

    xcorr_f = fft(sig1,Nfft).*conj(fft(sig2,Nfft));
    phase_f = unwrap(angle(fftshift(xcorr_f)), phase_wrap_step);
    H = [fkeep_freqs.' ones(length(fkeep_freqs),1)]; % used for LS estimation of delay
    p = phase_f(occbw_idx);
    m = H\p;
    ls_phase_est = H*m - m(2);
    ls_mean_phase_est = m(2);
    delay_phase_est = -Nfft*ls_phase_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    mean_delay_est = mean(delay_phase_est);
end

function plot3d_slice(f_vec, alpha_vec, z_mat, sym_rate)
    [Nalpha, Nfreq] = size(z_mat);
    alpha_max = max(alpha_vec);
    alpha_tol = 0.01;
    plot_opts = {'r-','g-','m-','y-'};
    
    figure
    hold on
    kk = 1;
    for i = 1:Nalpha
        if nargin == 4
            if(i == 1)
                plot_opt = 'b-';
                lw = 1.5;
            elseif (mod(alpha_vec(i),sym_rate) < alpha_tol*sym_rate)
                plot_opt = plot_opts{kk};
                lw = 1.5;
                kk = kk + 1;
            else
                plot_opt = 'k-';
                lw = 0.5;
            end
            plot3(f_vec, alpha_vec(i)*ones(Nfreq,1), z_mat(i,:), plot_opt, 'linewidth', lw)
            plot3([-sym_rate/2 -sym_rate/2], [0 alpha_max], [0 0], 'k--') % highlights occupied bandwidth in xy plane
            plot3([sym_rate/2 sym_rate/2], [0 alpha_max], [0 0], 'k--')
            
            plot3([-sym_rate/2 -sym_rate/2], [0 0], zlim, 'k--') % highlights occupied bandwidth in xz plane
            plot3([sym_rate/2 sym_rate/2], [0 0], zlim, 'k--')
            
            plot3(xlim, [sym_rate sym_rate], [0 0], 'k--') % highlights first cycle freq
            plot3([0 0], [0 0], zlim, 'k-') % highlights spectral freq 0 xy plane
            plot3([0 0], ylim, [0 0], 'k-') % highlights spectral freq 0 xz plane
        else
            plot3(f_vec, alpha_vec(i)*ones(Nfreq,1), z_mat(i,:))
        end
        
    end
    view(3)
    set(gca,'Ydir','reverse')
    grid on
    xlabel('f (spectral freq)', 'Rotation',+15)
    ylabel('\alpha (cycle freq)', 'Rotation',-25)
    zlabel('|SCD|')
    xticks(-0.5:0.2:0.5)
end
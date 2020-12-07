clear; close all
%% Simulation parameters
N = 100;     % length of first sequence (symbols)
sps = 4;   % samples per symbol for shaping filter
fbw = 7e6;  % occupied bandwidth of transmitter
% fs = fbw*sps; % 14e6;  % receiver sampling rate (Hz)
fs = 28e6;  % receiver sampling rate (Hz)
Ts = 1/fs;  % sample period (s)
D = 3*(2*rand-1)*Ts;  % delay (s)
sig = 0.0005; % standard deviation of noise
multi_sps = 3;  % oversample factor when multipath is added
M = 2;      % modulation order
pulse_shaping = 0;  % pulse shape using root raised cosine 0/1 - no/yes
span = 5;   % length of pulse shaping filter (symbols)
beta = 0.18; % excess bandwidth of pulse shaping filter

% Cyclostationary plot parameters (normalized units)
os_fac    = 2; % oversampling amount
os_fac_conj = 2;
sym_rate = fbw/fs; % noramalized symbol rate
alpha_scd = 0.:sym_rate/4:1; % define cycle frequencies to compute
alpa_sym_idx = find(alpha_scd == sym_rate);
% alpha_scd = [0 0.3 0.5 0.7];
freq_res = 0.15;
wind_o_n  = os_fac*round(freq_res*N*sps); % number of window taps
wind_o    = ones(wind_o_n,1)/wind_o_n; % normalized window

%% generate the signals
[tx_sig, raw_symbols, shaped_samples] = ...
    get_tx_signal(N, M, sps, fs, fbw, pulse_shaping, span, beta); 

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
plot(real(raw_symbols), 'b.-'); hold on
plot(imag(raw_symbols), 'r.-')
xlabel('Sample Number')
ylabel('Amplitude')
title('Input')
axis([0 inf -inf inf])

subplot(4,1,2)
plot(real(shaped_samples), 'b.-')
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
plot(faxis*1e-6, tx_sig_f, '.-'); hold on
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlabel('Frequency (MHz)')
ylabel('|fft|^2')
title('Signal Power Spectrum')

%% Cyclostationarity check of input sequence
nfft      = 2^(nextpow2(os_fac*length(shaped_samples)));
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
pdm = 0.15; % peak distance multipier just for plot display, what percent 
            % of samples around the main peak do you want to display
xtick_step = 16*Ts*1e9;
num_rows = multi_sps+2;

[z, lags] = xcorr(tx_sig,tx_sig);
[pval, pidx] = max(abs(z));
z_skew = skewness(abs(z(pidx-sps:pidx+sps)));
Nfft = max(1024, 2^(nextpow2(length(z)))); % must be even for correct results
faxis = fs*(-0.5:1/Nfft:0.5-1/Nfft);
tx_sig_f = 10*log10(fftshift(1/Nfft*abs(fft(tx_sig,Nfft)).^2));
rx_sig_no_multi_f = 20*log10(fftshift(1/Nfft*abs(fft(rx_sig_no_multi.',Nfft)).^2));
rx_sig_frac_multi_f = 20*log10(fftshift(1/Nfft*abs(fft(rx_sig_frac_multi.',Nfft)).^2));
rx_sig_int_multi_f = 20*log10(fftshift(1/Nfft*abs(fft(rx_sig_int_multi.',Nfft)).^2));

figure
subplot(num_rows,2,1)
plot_xcorr_time(z, N, pval, pidx, lags, Ts, sps, pdm, xtick_step)
xlims = xlim;
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
title(sprintf('Fractional Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))
  
subplot(num_rows,2,2)
% estimate the delay using phase
bw = 0.9*fbw; %fs/sps2; % the bandwidth to use for computations
removal_band = 0.2e6; % band around zero to ignore, no delay information at k=0
phase_wrap_step = pi;
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(tx_sig, tx_sig, faxis, Nfft, fs, bw, removal_band, phase_wrap_step);

plot_phase_frequency(phase_f, array_idx, faxis, fbw, m, p_est, mean_delay_est, 0)

subplot(num_rows,2,3)
[z, lags] = xcorr(rx_sig_no_multi,tx_sig);
[pval, pidx] = max(abs(z));
z_skew = skewness(abs(z(pidx-sps:pidx+sps)));

plot_xcorr_time(z, N, pval, pidx, lags, Ts, sps, pdm, xtick_step)
xlims = xlim;
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

subplot(num_rows,2,4)
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_no_multi, tx_sig, faxis, Nfft, fs, bw, removal_band, phase_wrap_step);

plot_phase_frequency(phase_f, array_idx, faxis, fbw, m, p_est, mean_delay_est, 1)

kk = 1;
for ii = 1:2:2*(num_rows-2)
    % plot correlation in time
    subplot(num_rows,2,ii+4)
    [z, lags] = xcorr(rx_sig_frac_multi(kk,:).',tx_sig);
    [pval, pidx] = max(abs(z));
    z_skew = skewness(abs(z(pidx-sps:pidx+sps)));
    
    plot_xcorr_time(z, N, pval, pidx, lags, Ts, sps, pdm, xtick_step)
    xlims = xlim;
    text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
    title(sprintf('Second Path Delay: %3.1f ns', kk*Ts/multi_sps*1e9))
    
    % plot phase of correlation in frequency
    subplot(num_rows,2,ii+5)
    [mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_frac_multi(kk,:).', tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);

    plot_phase_frequency(phase_f, array_idx, faxis, fbw, m, p_est, mean_delay_est, 1)
        
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

plot_mag_phase_together(tx_sig_f, phase_f-m, faxis, fbw, 0)
title(sprintf('Fractional Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))

subplot(num_rows,1,2)
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_no_multi, tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);    

plot_mag_phase_together(rx_sig_no_multi_f, phase_f-m, faxis, fbw, 1)
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

kk = 1;
for ii = 1:1:1*(num_rows-2)
    subplot(num_rows,1,ii+2)   
    [mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_frac_multi(kk,:).', tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);   

    plot_mag_phase_together(rx_sig_frac_multi_f(:,ii), phase_f, faxis, fbw, 1)
    title(sprintf('Second Path Delay: %3.1f ns', kk*Ts/multi_sps*1e9))
        
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

pdm = 0.15; % peak distance multipier
num_rows = ceil(int_max_delay/step1)+2;

figure
subplot(num_rows,2,1)
plot_xcorr_time(z, N, pval, pidx, lags, Ts, sps, pdm, xtick_step)
xlims = xlim;
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
title(sprintf('Integer Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))
  
subplot(num_rows,2,2)
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(tx_sig, tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);  
plot_phase_frequency(phase_f, array_idx, faxis, fbw, m, p_est, mean_delay_est, 1)

subplot(num_rows,2,3)
[z, lags] = xcorr(rx_sig_no_multi,tx_sig);
[pval, pidx] = max(abs(z));
z_skew = skewness(abs(z(pidx-sps:pidx+sps)));
plot_xcorr_time(z, N, pval, pidx, lags, Ts, sps, pdm, xtick_step)
xlims = xlim;
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

subplot(num_rows,2,4)
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_no_multi, tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);  
plot_phase_frequency(phase_f, array_idx, faxis, fbw, m, p_est, mean_delay_est, 1)

kk = 1;
for ii = 1:2:2*(num_rows-2)
    subplot(num_rows,2,ii+4)
    [z, lags] = xcorr(rx_sig_int_multi(kk,:).',tx_sig);
    [pval, pidx] = max(abs(z));
    z_skew = skewness(abs(z(pidx-sps:pidx+sps)));
    plot_xcorr_time(z, N, pval, pidx, lags, Ts, sps, pdm, xtick_step)
    xlims = xlim;
    text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
    title(sprintf('Second Path Delay: %3.1f ns', delays(kk)*Ts*1e9))
    
    subplot(num_rows,2,ii+5)
    [mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_int_multi(kk,:).', tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);  
    plot_phase_frequency(phase_f, array_idx, faxis, fbw, m, p_est, mean_delay_est, 1)
        
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

plot_mag_phase_together(tx_sig_f, phase_f, faxis, fbw, 0)
title(sprintf('Integer Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))

subplot(num_rows,1,2)
[mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
    get_phase_est(rx_sig_no_multi, tx_sig, faxis, Nfft, fs, ...
    bw, removal_band, phase_wrap_step);  

plot_mag_phase_together(rx_sig_no_multi_f, phase_f, faxis, fbw, 1)
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

kk = 1;
for ii = 1:1:1*(num_rows-2)
    subplot(num_rows,1,ii+2)    
    [mean_delay_est, phase_f, p_est, m, fk, array_idx, kidx] = ...
        get_phase_est(rx_sig_int_multi(kk,:).', tx_sig, faxis, Nfft, fs, ...
        bw, removal_band, phase_wrap_step); 
    
    plot_mag_phase_together(rx_sig_int_multi_f(:,ii), phase_f, faxis, fbw, 1)
    title(sprintf('Second Path Delay: %3.1f ns', delays(kk)*Ts*1e9))
    
    kk = kk + 1;
end

%% 
%%%%%%%%%%%%%%%%%%%% FUNCTION DEFS BELOW HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tx_sig, raw_symbols, shaped_samples] = ...
    get_tx_signal(N, M, sps, fs, fbw, pulse_shaping, span, beta)

    % create the baseline modulated signal
    xb = randi([0 M-1], log2(M)*N, 1);
    raw_symbols = qammod(xb, M, 'gray'); % raw symbols
    x2 = upsample(raw_symbols, sps);
    if pulse_shaping == 1
        rrc = rcosdesign(beta, span, sps, 'sqrt');
        shape_filt = rrc.'/max(rrc);
    else
        shape_filt = ones(sps,1);
    end
    shaped_samples = filter(shape_filt, 1, x2); % shaped samples

    [P,Q] = rat((fs/fbw)/sps);
    if P*Q == 1
        tx_sig = shaped_samples/sqrt(mean(abs(shaped_samples).^2)); % normalize power to 1
    else
        x4 = resample(shaped_samples, P, Q);
        tx_sig = x4/sqrt(mean(abs(x4).^2)); % normalize power to 1
    end
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
    kidx = faxis*Nfft/fs;
    fkeep_idx = find(faxis > -bw/2 & faxis < bw/2);
    frem_idx = find(faxis > -removal_band & faxis < removal_band);
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
            % use a different color and  linewidth for each cycle freq
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

function plot_xcorr_time(z, N, pval, pidx, lags, Ts, sps, pdm, xtick_step)
    plot(lags*Ts*1e9, abs(z), '.-'); hold on
    plot(lags(pidx)*Ts*1e9, pval, 'rx')
    plot(lags(pidx-sps:pidx+sps)*Ts*1e9, abs(z(pidx-sps:pidx+sps)), 'o', 'MarkerSize', 4)
    axis([-pdm*N*sps*Ts*1e9 pdm*N*sps*Ts*1e9 0 1.2*max(abs(z))])
    xlims = xlim;
    xticks([xlims(1):xtick_step:xlims(2)])
%     text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
    xlabel('Lag (ns)')
    ylabel('|xcorr|')
%     title(sprintf('Fractional Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))
end

function plot_phase_frequency(phase_f, array_idx, faxis, fbw, m, p_est, ...
    mean_delay_est, tight)

    plot(faxis*1e-6, phase_f-m, 'b-'); hold on
    plot(faxis(array_idx)*1e-6, phase_f(array_idx)-m, 'b.')
    plot(faxis(array_idx)*1e-6, p_est, 'b-')
    if tight == 1
        axis tight
    else
        axis([faxis(1)*1e-6 faxis(end)*1e-6 -1 1])
    end
    plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
    plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
    xlims = xlim;
    ylims = ylim;
    text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f', mean_delay_est))
    ylabel('Phase')
    xlabel('Frequency (MHz)')
end

function plot_mag_phase_together(mag_f, phase_f, faxis, fbw, tight)
    yyaxis left
    plot(faxis*1e-6, mag_f, 'b-'); hold on
    axis([-inf inf -100 0])
    ylabel('PSD')
    yyaxis right
    plot(faxis*1e-6, phase_f, 'r-')
    if tight == 1
        axis tight
    else
        axis([faxis(1)*1e-6 faxis(end)*1e-6 -1 1])
    end
    plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
    plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
    xlims = xlim;
    ylims = ylim;
    ylabel('Phase')
    xlabel('Frequency (MHz)')
%     title(sprintf('Fractional Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))
end
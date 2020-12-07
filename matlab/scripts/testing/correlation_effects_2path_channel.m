clear; close all

N = 30;     % length of first sequence (symbols)
fs = 10e6;  % receiver sampling rate (Hz)
fbw = 7e6;  % occupied bandwidth of transmitter
Ts = 1/fs;  % sample period (s)
D = 3*(2*rand-1)*Ts;  % delay (s)
sig = 0.05; % standard deviation of noise
sps1 = 3;  % samples per sample at high rate when multipath is added
sps2 = 4;   % samples per symbol for shaping filter
M = 4;      % modulation order
span = 5;   % length of shaping filter (symbols)
beta = 0.2; % excess bandwidth of shaping filter

% Cyclostationary plot parameters
alpha_scd = 0.:0.1:0.7; % define cycle frequencies to compute
alpha_scd = [0 0.3 0.5 0.7];
wind_o_n  = 20; % number of window taps
wind_o    = ones(wind_o_n,1)/wind_o_n; % normalized window
os_fac    = 10; % oversampling amount
os_fac_conj = 10;

% create the baseline modulated signal
xb = randi([0 M-1], log2(M)*N, 1);
x1 = qammod(xb, M, 'gray');
x2 = upsample(x1, sps2);
rrc = rcosdesign(beta, span, sps2, 'sqrt');
rrc = rrc.'/max(rrc);
shape_filt = ones(sps2,1);
x3 = conv(rrc, x2);
% x3 = conv(shape_filt, x2);

[P,Q] = rat((fs/fbw)/sps2);
x4 = resample(x3, P, Q);
% x4=x3;
% x4 = sin(2*pi*1/sps2*(0:sps2*N-1)).';
x5 = x4/sqrt(mean(abs(x4).^2)); % normalize power to 1

% delay the sequence
max_num_delay_samps = ceil(D/Ts)-1;
y1 = delayseq([x5; zeros(max_num_delay_samps,1)], D, fs);
y1 = y1 + sig*(randn(length(y1),1)+1j*randn(length(y1),1));

% add fractional delayed multipaths
y5 = zeros(sps1, length(y1)+2);
for ii = 1:sps1
    y2 = resample(y1, sps1, 1);
    h = [1, zeros(1, ii-1), 0.95];
    hc{ii} = h;
    y3 = conv(y2,h);
    y4 = resample(y3,1,sps1);
    y5(ii,1:length(y4)) = y4;
end

% add integer delayed multipaths
int_max_delay = 11; % samples
step1 = 4;
y7 = zeros(ceil(int_max_delay/step1), length(y1)+int_max_delay);
kk = 1;
for ii = 1:step1:int_max_delay
    delays(kk) = ii;
    h = [1 zeros(1, ii-1) 0.95];
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
% plot(real(x5), 'b.-'); hold on
plot(real(x3), 'b.-')
xlabel('Sample Number')
ylabel('Amplitude')
title('Shaped Channel Input')
axis([0 inf -inf inf])

subplot(4,1,3)
plot(real(y1), 'b.-')
xlabel('Sample Number')
ylabel('Amplitude')
title('Delayed AWGN Channel Output')
axis([0 inf -inf inf])

subplot(4,1,4)
Nfft = max(1024, 2^(nextpow2(size(y7,2)))); % must be even for correct results
faxis = fs*(-0.5:1/Nfft:0.5-1/Nfft);
x5_f = 10*log10(fftshift(1/1024*abs(fft(x5,Nfft)).^2));
y1_f = 20*log10(fftshift(1/1024*abs(fft(y1.',Nfft)).^2));
y5_f = 20*log10(fftshift(1/1024*abs(fft(y5.',Nfft)).^2));
y7_f = 20*log10(fftshift(1/1024*abs(fft(y7.',Nfft)).^2));
lh = [];
lh(1) = plot(faxis*1e-6, x5_f, '.-'); hold on
% lht = plot(faxis*1e-6, y5_f);
lht = plot(faxis*1e-6, y7_f);
lh(2) = lht(1);
plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
xlabel('Frequency (MHz)')
ylabel('|fft|^2')
title('Signal Power Spectrum')
legend(lh, 'Transmitted','Received')

%% Cyclostationarity check of input sequence
nfft      = os_fac*length(x3);
scd_fsm_mat = os_fac*cyclo_scd_fsm(x3,alpha_scd,wind_o,nfft);
scd_conj_fsm_mat = os_fac_conj*cyclo_scd_fsm_conj(x3,alpha_scd,wind_o,nfft);

norm_f_vec = (1:nfft) - nfft/2;
norm_f_vec = norm_f_vec/max(2*abs(norm_f_vec));
figure
for ind = 1:size(scd_fsm_mat,1)
   lh(ind) = plot(norm_f_vec,10*log10(abs(scd_fsm_mat(ind,:)))); hold on 
%    plot(faxis/fs, x5_f, '.-')
end
hold off
xlabel('Normalized Frequency','FontSize',14);
ylabel('|SCD|','FontSize',14);
grid minor;
grid on;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
axis tight;
ylim([-20 10]);
xlim([-0.4 0.4])
lgd_han = legend(lh, string(alpha_scd));
title("Cyclic Frequency (\alpha)");
title("Frequency Smoothed SCD");

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
phase_wrap_step = pi;
xtick_step = 4*Ts*1e9;
num_rows = sps1+2;
figure
subplot(num_rows,2,1)
plot(lags*Ts*1e9, abs(z), '.-'); hold on
plot(lags(pidx)*Ts*1e9, pval, 'rx')
plot(lags(pidx-sps2:pidx+sps2)*Ts*1e9, abs(z(pidx-sps2:pidx+sps2)), 'o', 'MarkerSize', 4)
axis([-pdm*N*sps2*Ts*1e9 pdm*N*sps2*Ts*1e9 0 1.2*max(abs(z))])
xlims = xlim;
xticks([xlims(1):xtick_step:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Fractional Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))
  
subplot(num_rows,2,2)
z_f = fft(x5,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);

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
xticks([xlims(1):xtick_step:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

subplot(num_rows,2,4)
% z_f = fft(ifftshift(z));
z_f = fft(y1,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);

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
    xticks([xlims(1):xtick_step:xlims(2)])
    text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
    text(xlims(1), 0.75*pval, sprintf('Channel: h = [%s]', num2str(hc{kk},'%2.1f ')))
    xlabel('Lag (ns)')
    ylabel('|xcorr|')
    title(sprintf('Second Path Delay: %3.1f ns', kk*Ts/sps1*1e9))
    
%     z_f = fft(ifftshift(z));
    z_f = fft(y5(kk,:).',Nfft).*conj(fft(x5,Nfft));
    phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);
    
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

%% Fractional Delay Magnitude ontop of Phase Plots
pdm = 0.15; % peak distance multipier for plots
phase_wrap_step = pi;
xtick_step = 4*Ts*1e9;
num_rows = sps1+2;

figure
subplot(num_rows,1,1)
  
z_f = fft(x5,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), pi);

p = phase_f(array_idx);
m = H\p;
p_est = H*m;
delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
mean_delay_est = mean(delay_phase_est);

delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
yyaxis left
plot(faxis*1e-6, x5_f(:,1), 'b-'); hold on
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
z_f = fft(y1,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), pi);

p = phase_f(array_idx);
m = H\p;
p_est = H*m-m(2);
delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
mean_delay_est = mean(delay_phase_est);

delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
yyaxis left
plot(faxis*1e-6, y1_f, 'b-'); hold on
axis([-inf inf -100 0])
ylabel('PSD')
yyaxis right
plot(faxis*1e-6, phase_f-m(2), 'r-')
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
    z_f = fft(y5(kk,:).',Nfft).*conj(fft(x5,Nfft));
    phase_f = unwrap(angle(fftshift(z_f)), pi);
    
    p = phase_f(array_idx);
    m = H\p;
    p_est = H*m-m(2);
    delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    mean_delay_est = mean(delay_phase_est);
    
    delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    yyaxis left
    plot(faxis*1e-6, y5_f(:,ii), 'b-'); hold on
    axis([-inf inf -100 0])
    ylabel('PSD')
    yyaxis right
    plot(faxis*1e-6, phase_f-m(2), 'r-')
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
% SCD or conj SCD. Use os_fac 2 for conj SCD.
scd_fsm_mat = os_fac*cyclo_scd_fsm(x5,alpha_scd,wind_o,nfft);
scd_conj_fsm_mat = os_fac_conj*cyclo_scd_fsm_conj(x5,alpha_scd,wind_o,nfft);

norm_f_vec = (1:nfft) - nfft/2;
norm_f_vec = norm_f_vec/max(2*abs(norm_f_vec));

% figure
% surf(norm_f_vec, alpha_scd, abs(scd_fsm_mat))
% t = 0:pi/50:10*pi;
% st = sin(t);
% ct = cos(t);
% plot3(st,ct,t,'.')

figure
subplot(num_rows,2,1)
for ind = 1:size(scd_fsm_mat,1)
   lh(ind) = plot(norm_f_vec,10*log10(abs(scd_fsm_mat(ind,:)))); hold on 
%    plot(faxis/fs, x5_f, '.-')
end
hold off
xlabel('Normalized Frequency','FontSize',14);
ylabel('|SCD|','FontSize',14);
grid minor;
grid on;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
axis tight;
ylim([-20 10]);
xlim([-0.4 0.4])
lgd_han = legend(lh, string(alpha_scd));
title(lgd_han, "Cyclic Frequency (\alpha)");
title("Frequency Smoothed SCD Tx");

subplot(num_rows,2,2)
for ind = 1:size(scd_conj_fsm_mat,1)
   plot(norm_f_vec,10*log10(abs(scd_conj_fsm_mat(ind,:)))); 
   hold on 
end
hold off
xlabel('Normalized Frequency','FontSize',14);
ylabel('|SCD|','FontSize',14);
grid minor;
grid on;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
axis tight;
ylim([-20 10]);
xlim([-0.4 0.4])
lgd_han = legend(string(alpha_scd));
title(lgd_han,"Cyclic Frequency (\alpha)");
title("Frequency Smoothed Conjugate SCD");

nfft      = os_fac*length(y1);
scd_fsm_mat = os_fac*cyclo_scd_fsm(y1,alpha_scd,wind_o,nfft);
scd_conj_fsm_mat = os_fac_conj*cyclo_scd_fsm_conj(y1,alpha_scd,wind_o,nfft);

norm_f_vec = (1:nfft) - nfft/2;
norm_f_vec = norm_f_vec/max(2*abs(norm_f_vec));

subplot(num_rows,2,3)
norm_f_vec = (1:nfft) - nfft/2;
norm_f_vec = norm_f_vec/max(2*abs(norm_f_vec));
for ind = 1:size(scd_fsm_mat,1)
   plot(norm_f_vec,10*log10(abs(scd_fsm_mat(ind,:)))); 
   hold on 
end
hold off
xlabel('Normalized Frequency','FontSize',14);
ylabel('|SCD|','FontSize',14);
grid minor;
grid on;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
axis tight;
ylim([-20 10]);
xlim([-0.4 0.4])
% lgd_han = legend(string(alpha_scd));
% title(lgd_han,"Cyclic Frequency (\alpha)");
title("Frequency Smoothed SCD y1");

subplot(num_rows,2,4)
for ind = 1:size(scd_conj_fsm_mat,1)
   plot(norm_f_vec,10*log10(abs(scd_conj_fsm_mat(ind,:)))); 
   hold on 
end
hold off
xlabel('Normalized Frequency','FontSize',14);
ylabel('|SCD|','FontSize',14);
grid minor;
grid on;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
axis tight;
ylim([-20 10]);
xlim([-0.4 0.4])
% lgd_han = legend(string(alpha_scd));
% title(lgd_han,"Cyclic Frequency (\alpha)");
title("Frequency Smoothed Conjugate SCD");

kk = 1;
for ii = 1:2:2*(num_rows-2)
    nfft      = os_fac*length(y5(kk,:).');
    scd_fsm_mat = os_fac*cyclo_scd_fsm(y5(kk,:).',alpha_scd,wind_o,nfft);
    
    norm_f_vec = (1:nfft) - nfft/2;
    norm_f_vec = norm_f_vec/max(2*abs(norm_f_vec));
    
    subplot(num_rows,2,ii+4)
    for ind = 1:size(scd_fsm_mat,1)
       plot(norm_f_vec,10*log10(abs(scd_fsm_mat(ind,:)))); 
       hold on 
    end
    hold off
    xlabel('Normalized Frequency','FontSize',14);
    ylabel('|SCD|','FontSize',14);
    grid minor;
    grid on;
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    axis tight;
    ylim([-20 10]);
    xlim([-0.4 0.4])
%     lgd_han = legend(string(alpha_scd));
%     title(lgd_han,"Cyclic Frequency (\alpha)");
    title(sprintf("Frequency Smoothed SCD 2nd Path Delay %3.0f ns",kk*Ts/sps1*1e9));
    
    scd_conj_fsm_mat = os_fac_conj*cyclo_scd_fsm_conj(y5(kk,:).',alpha_scd,wind_o,nfft);
    subplot(num_rows,2,ii+5)
    for ind = 1:size(scd_conj_fsm_mat,1)
       plot(norm_f_vec,10*log10(abs(scd_conj_fsm_mat(ind,:)))); 
       hold on 
    end
    hold off
    xlabel('Normalized Frequency','FontSize',14);
    ylabel('|SCD|','FontSize',14);
    grid minor;
    grid on;
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    axis tight;
    ylim([-20 10]);
    xlim([-0.4 0.4])
%     lgd_han = legend(string(alpha_scd));
%     title(lgd_han,"Cyclic Frequency (\alpha)");
    title("Frequency Smoothed Conjugate SCD");
        
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
xticks([xlims(1):xtick_step:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Integer Delayed Multipaths, Ts = %3.0f ns\nAutocorrelation', Ts*1e9))
  
subplot(num_rows,2,2)
z_f = fft(x5,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);

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
text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f ns', mean_delay_est))
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
xticks([xlims(1):xtick_step:xlims(2)])
text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
xlabel('Lag (ns)')
ylabel('|xcorr|')
title(sprintf('Direct Path Delay: %3.1f ns', D*1e9))

subplot(num_rows,2,4)
% z_f = fft(ifftshift(z));
z_f = fft(y1,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);

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
text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f ns', mean_delay_est))
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
    xticks([xlims(1):xtick_step:xlims(2)])
    text(xlims(1), pval, sprintf('Skewness: %4.2f', z_skew))
    text(xlims(1), 0.75*pval, sprintf('Channel: h = [%s]', num2str(hc{kk},'%2.1f ')))
    xlabel('Lag (ns)')
    ylabel('|xcorr|')
    title(sprintf('Second Path Delay: %3.1f ns', delays(kk)*Ts*1e9))
    
%     z_f = fft(ifftshift(z));
    z_f = fft(y7(kk,:).',Nfft).*conj(fft(x5,Nfft));
    phase_f = unwrap(angle(fftshift(z_f)), phase_wrap_step);
    
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
    text(0, ylims(2)*0.7, sprintf('Est Delay: %4.2f ns', mean_delay_est))
    ylabel('Phase')
    xlabel('Frequency (MHz)')
        
    kk = kk + 1;
end

%% Integer Delay Magnitude ontop of Phase Plots
pdm = 0.15; % peak distance multipier for plots
phase_wrap_step = pi;
xtick_step = 4*Ts*1e9;
num_rows = sps1+2;

figure
subplot(num_rows,1,1)
  
z_f = fft(x5,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), pi);

p = phase_f(array_idx);
m = H\p;
p_est = H*m;
delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
mean_delay_est = mean(delay_phase_est);

delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
yyaxis left
plot(faxis*1e-6, x5_f(:,1), 'b-'); hold on
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
z_f = fft(y1,Nfft).*conj(fft(x5,Nfft));
phase_f = unwrap(angle(fftshift(z_f)), pi);

p = phase_f(array_idx);
m = H\p;
p_est = H*m-m(2);
delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
mean_delay_est = mean(delay_phase_est);

delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
yyaxis left
plot(faxis*1e-6, y1_f, 'b-'); hold on
axis([-inf inf -100 0])
ylabel('PSD')
yyaxis right
plot(faxis*1e-6, phase_f-m(2), 'r-')
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
    z_f = fft(y7(kk,:).',Nfft).*conj(fft(x5,Nfft));
    phase_f = unwrap(angle(fftshift(z_f)), pi);
    
    p = phase_f(array_idx);
    m = H\p;
    p_est = H*m-m(2);
    delay_phase_est = -Nfft*p_est./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    mean_delay_est = mean(delay_phase_est);
    
    delay_phase = -Nfft*phase_f(array_idx)./(2*pi*fs*kidx.')*1e9; % estimated delay (ns)
    yyaxis left
    plot(faxis*1e-6, y7_f(:,ii), 'b-'); hold on
    axis([-inf inf -100 0])
    ylabel('PSD')
    yyaxis right
    plot(faxis*1e-6, phase_f-m(2), 'r-')
    axis tight
    plot([-fbw/2 -fbw/2]*1e-6, ylim, 'k--')
    plot([fbw/2 fbw/2]*1e-6, ylim, 'k--')
    xlims = xlim;
    ylims = ylim;
    ylabel('Phase')
    xlabel('Frequency (MHz)')
%     title(sprintf('Channel %i, Ntap %i, Npath %i, MaxAmp %.1f', kk, Ntap, Npath, max_nlos_amp))
%     title(sprintf('Channel %3.0f, Indices %s', kk, string(hidx_sv{kk})))
    
    kk = kk + 1;
end
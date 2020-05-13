clear; close all

%% Sim parameters
perturb_est = 1;    % use imperfect estimates of delay and attenuation
frac_delay = 1; % use fractional delays or integer multiples of sample rate
Ntrials = 30;   % number of trials for statistical results
Nsym = 100;     % number of symbols making up tx signal
Nrx = 3;        % number of receivers
Nmp = 2;        % number of multipaths at all but one receiver
snrdb = 20;     % SNR at each receiver
sps = 2;        % samples per symbol at tx
frx = 100e6;    % receiver sample rate (Hz)
delay_spread = 500e-9;   % time difference between first received path and last (s)
beta = 0.4;     % tx pulse shaping excess bandwidth
span = 10;      % tx pulse shaping filter span in symbols
c = 299792458;  % speed of light (m/s)

%% Generate signal
x1 = 2*randi([0 1], Nsym, 1)-1 + 1j*(2*randi([0 1], Nsym, 1)-1); % qpsk symbols
x2 = upsample(x1,sps);
rrc = rcosdesign(beta, span, sps); % pulse shaping
rrc = rrc.'/max(rrc);
x3 = conv(rrc, x2);
x3 = x3/sqrt(mean(abs(x3).^2)); % normalize power to 1
Nsamps = length(x3);

%% Generate a set of time delays and complex attenuations
max_toa = 200;  % max toa (ns)
Ts = 1/frx;
max_toa_samps = ceil(max_toa*frx*1e-9);
if frac_delay == 0
    toas_idx = sort(randperm(max_toa_samps, Nrx-1),2);
    toas_los = toas_idx/frx; % los delay from first receiver
else
    toas_los = sort(max_toa*rand(1,Nrx-1),2)*1e-9;
end

%% Delay received signals relative to the first received signal
x4 = zeros(max_toa_samps+Nsamps,Nrx);
x4(1:Nsamps,1) = x3; % first recieiver signal is not delayed
for ii = 2:Nrx
    xd = delayseq([x3; zeros(max_toa_samps,1)], toas_los(ii-1), frx);
    x4(:,ii) = [xd; zeros((max_toa_samps+Nsamps)-length(xd),1)];
end
% plot(real(x4))

%% Estimate the tdoa before adding multipath
% Cross correlate all receivers against the template
corr_out0 = zeros(2*size(x4,1)-1, Nrx-1);
lags0 = zeros(2*size(x4,1)-1, Nrx-1);
x4means = mean(x4);
x4 = x4 - x4means; % make it zero mean
for ii = 2:Nrx
    [corr_out0(:,ii-1), lags0(:,ii-1)] = xcorr(x4(:,ii), x4(:,1));
end
corr_mag_sq0 = abs(corr_out0);
[xmax0, idx_max0] = max(corr_mag_sq0);
toas_pre_est = lags0(idx_max0,1)*1/frx;

%% Estimate the fine tdoa after SRC
samps_from_peak = 1;
for ii = 2:Nrx
    w = lags0(idx_max0(ii-1)-samps_from_peak:idx_max0(ii-1)+samps_from_peak,ii-1);
    z = corr_mag_sq0(idx_max0(ii-1)-samps_from_peak:idx_max0(ii-1)+samps_from_peak,ii-1);
%     plot(w,z,'.-') % debug

    % Parabolic Interpolation
    H = [w.^2 w ones(length(w),1)];
    theta = H\z; % the coefficients of the line fitting the data points

    wm = -theta(2)/(2*theta(1));        % value of x where y is max
    wm_t = wm*Ts*1e9;
    toas_pre_refined_est(ii-1) = wm_t*1e-9;
end

toas_los_coarse_est = zeros(Nrx-1,Ntrials);
toas_mplos_est = zeros(Nrx-1,Ntrials);
toas_los_est_refined = zeros(Nrx-1,Ntrials);
for nn = 1:Ntrials
    %% Add multipath
    num_taps = ceil(delay_spread*frx); % number of samples per delay spread interval
    for ii = 2:Nrx
        toas_multi(:,ii-1) = sort(randperm(num_taps, max(0,Nmp-1))); % num_paths-1 bc one path is the direct path
        attens_multi(:,ii-1) = rand(1,Nmp) .* exp(1j*2*pi*rand(1,Nmp));

        h = zeros(num_taps,1);
        h(1) = attens_multi(1,ii-1); % direct path uneffected here
        h(toas_multi(:,ii-1)) = attens_multi(2:end,ii-1);

        m1(:,ii) = conv(x4(:,ii),h);
    end
    m1(1:length(x4(:,1)),1) = x4(:,1);

    %% Add noise
    noise_var = 10^(-snrdb/10);
    w = sqrt(noise_var/2)*(randn(size(m1))+1j*randn(size(m1)));
    r1 = m1 + w;

    %% Estimate the tdoa before multipath removal
    % Cross correlate all receivers against the template
    corr_out = zeros(2*size(r1,1)-1, Nrx-1);
    lags1 = zeros(2*size(r1,1)-1, Nrx-1);
    r1means = mean(r1);
    r1 = r1 - r1means; % make it zero mean
    for ii = 2:Nrx
        [corr_out1(:,ii-1), lags1(:,ii-1)] = xcorr(r1(:,ii), r1(:,1));
    end
    corr_mag_sq1 = abs(corr_out1);
    [xmax1, idx_max1] = max(corr_mag_sq1);
    toas_mplos_est(:,nn) = lags1(idx_max1,1)*1/frx;

    %% Successive receiver cancellation (SRC) technique
    % 1 - Determine which receiver will act as the template
    t_idx_est = 1;  % genie delivered template estimate

    % 2a - Estimate the delay between the template first receiver and others
    delay_est = round(toas_los*frx);   % genie delivered estimates
    delay_std = 1;
    if perturb_est == 1
        delay_est = delay_est + round(delay_std*randn);
    end

    % 2b = Estimate the attenuation complex multipliers
    atten_est = attens_multi(1,:);   % genie delivered estimates
    atten_std = 0.2;
    if perturb_est == 1
        atten_est = atten_est + atten_std*randn(1,Nrx-1);
    end

    % 3 - Now remove the template from the other receivers and then remove the
    % isolated NLOS from the receivers leaving a clean version of the signal
    r2 = src(r1, t_idx_est, delay_est, atten_est);

    %% Estimate the coarse tdoa after SRC
    % Cross correlate all receivers against the template
    corr_out2 = [];%zeros(2*size(r2,1)-1, Nrx-1);
    lags2 = [];%zeros(2*size(r1,1)-1, Nrx-1);
    r2means = mean(r2);
    r2 = r2 - r2means; % make it zero mean
    for ii = 2:Nrx
        [corr_out2(:,ii-1), lags2(:,ii-1)] = xcorr(r2(:,ii), r2(:,1));
    end
    corr_mag_sq2 = abs(corr_out2);
    [xmax2, idx_max2] = max(corr_mag_sq2);
    toas_los_coarse_est(:,nn) = lags2(idx_max2,1)*1/frx;
    
    %% Estimate the fine tdoa after SRC
    samps_from_peak = 1;
    for ii = 2:Nrx
        w = lags2(idx_max2(ii-1)-samps_from_peak:idx_max2(ii-1)+samps_from_peak,ii-1);
        z = corr_mag_sq2(idx_max2(ii-1)-samps_from_peak:idx_max2(ii-1)+samps_from_peak,ii-1);
%         plot(w,z,'.-') % debug

        % Parabolic Interpolation
        H = [w.^2 w ones(length(w),1)];
        theta = H\z; % the coefficients of the line fitting the data points

        wm = -theta(2)/(2*theta(1));        % value of x where y is max
        wm_t = wm*Ts*1e9;
        toas_los_refined_est(ii-1,nn) = wm_t*1e-9;
    end
end

%% Plots
figure
plot(real(x3)); hold all
plot(real(x4(:,2)))
title('Time Domain Delayed Signals')
xlabel('Sample Number')
ylabel('Amplitude')

figure
subplot(Nrx,2,1)
plot(real(r1(:,1))); hold all
plot(imag(r1(:,1)))
title('Template Rx1')

subplot(Nrx,2,3)
plot(real(r1(:,2))); hold all
plot(imag(r1(:,2)))
title('Multipath Corrupted Rx2')

subplot(Nrx,2,5)
plot(real(r1(:,3))); hold all
plot(imag(r1(:,3)))
title('Multipath Corrupted Rx3')

subplot(Nrx,2,2)
plot(real(r2(:,1))); hold all
plot(imag(r2(:,1)))
title('Template Unchanged Rx1')

subplot(Nrx,2,4)
plot(real(r2(:,2))); hold all
plot(imag(r2(:,2)))
title('Multipath Removed Rx2')

subplot(Nrx,2,6)
plot(real(r2(:,3))); hold all
plot(imag(r2(:,3)))
title('Multipath Removed Rx3')

figure
xmin = inf;
xmax = -inf;
for nn = 1:Nrx-1
    toas_los_coarse_esti = toas_los_coarse_est(nn,:).'*1e9;
    toas_los_refined_esti = toas_los_refined_est(nn,:).'*1e9;
    toas_mplos_esti = toas_mplos_est(nn,:).'*1e9;
    mean_tdoa = mean(toas_los_refined_esti);
    std_tdoa = std(toas_los_refined_esti);
    htrue = plot(toas_los(nn)*1e9, nn, 'kx', 'markersize', 8); hold all
    hmp = plot(toas_mplos_esti, nn, 'rd');
    hcoarse = plot(toas_los_coarse_esti, nn, 'ro');
    hrefine = plot(toas_los_refined_esti, nn, 'b.');
    xmin = min([xmin; toas_los(nn)*1e9; toas_los_coarse_esti; toas_los_refined_esti; toas_mplos_esti]);
    xmax = max([xmax; toas_los(nn)*1e9; toas_los_coarse_esti; toas_los_refined_esti; toas_mplos_esti]);
    axis([min(-1, 1.1*xmin) max(1.1*xmax,1) 0.5 Nrx])
    text(toas_los(nn)*1e9, nn, sprintf('\\mu = %3.2f, \\sigma = %3.2f', mean_tdoa, std_tdoa), ...
        'horizontalalignment', 'center', 'verticalalignment', 'top',...
        'FontSize', 8);
end
y_values = 1:Nrx-1;
ylabels = {'1,2', '1,3'};
set(gca, 'Ytick', y_values, 'YTickLabel',ylabels);
xlabel('TDOA (ns)')
ylabel('Receiver Correlation Pairs')
hleglines = [htrue(1) hmp(1) hcoarse(1) hrefine(1)];
legend(hleglines, 'True TDOA', 'No SRC TDOA Estimate', ...
    'SRC Coarse TDOA Estimate', 'SRC Fine TDOA Estimate', 'Location', 'North')
set(gcf, 'Position',  [100, 100, 1300, 800])
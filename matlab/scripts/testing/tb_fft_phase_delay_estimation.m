% this script generates a sampled tx signal that is then delayed at each 
% receiver. It plots the time domain and correlation of the receiver pairs
% to estimate the delay. 
clear; close all

%% Simulation parameters
Ts = 10;            % receiver sampling period (ns)
delays = [0 5 27 100];  % delay to add (ns) at each receiver, length 
                    % corresponds to number of rxs, leave first delay as 0
snrdb = inf;        % signal to noise ratio of received signal in dB
Ntrials = 20;       % number of trials to test estimation accuracy
Nsym = 100;         % transmitted sequence length in symbols
sps = 2;            % samples per symbol at transmitter
span = 10;          % symbol span of shaping filter
beta = 0.5;         % shaping excess bandwidth factor
delay_spread = 300; % time between first received path and last (ns)
multi_option = 0;   % 0 = no multipath, 1 = 2 paths, 2 = random number paths
nlos_max_pw = 2;    % the maximum multiple the NLOS path power can be compared to the LOS power

%% Generate a test waveform, shaped PRN sequence
x1 = 2*randi([0 1], Nsym, 1)-1;
x2 = upsample(x1,sps);
rrc = rcosdesign(beta, span, sps);
rrc = rrc.'/max(rrc);
x3 = conv(rrc, x2);
x4 = x3/sqrt(mean(abs(x3).^2)); % normalize power to 1
sigman = 1/sqrt(10^(snrdb/10)); % noise standard deviation
Nsamp = length(x4);

%% Create a delayed version of the signal at each receiver
num_rx = length(delays);
fs = 1/Ts;
max_num_delay_samps = ceil(max(delays)/Ts)-1;
Nsamp = max_num_delay_samps + Nsamp;
y2 = zeros(Nsamp,num_rx);  % will contain the received sequences
for ii = 1:num_rx
    y1 = delayseq([x4; zeros(max_num_delay_samps,1)], delays(ii), fs);
    y2(:,ii) = [y1; zeros((Nsamp)-length(y1),1)];
end

%%
Nfft = 2*Nsamp-1;
% corrs_t = zeros(Ntrials, Nfft, num_rx-1);
for nn = 1:Ntrials
    
    %% add multipath
    if multi_option == 1 % 2 paths
        % max number of samples that can fit into the delay spread amount
        num_taps = ceil(delay_spread*fs);
        for ii = 1:num_rx
            rand_idx = randperm(num_taps,1); % random tap index
            channel_taps = zeros(num_taps,1);
            channel_taps(1) = 1; % direct path
            channel_taps(rand_idx) = sqrt(nlos_max_pw)*rand.*exp(1j*2*pi*rand); % indirect path
        
            y3(:,ii) = conv(y2(:,ii),channel_taps);
        end    
    elseif multi_option == 2 % random num paths 
        num_taps = ceil(delay_spread*fs);
        for ii = 1:num_rx
            num_paths = randi([0 num_taps]);
            rand_idxs = sort(randperm(num_taps, max([0 num_paths-1]))).';
            channel_taps = zeros(num_taps,1);
            channel_vals = sqrt(nlos_max_pw)*rand(num_paths-1,1) .* ...
                exp(1j*2*pi*rand(num_paths-1,1)); % indirect path
            channel_taps(1) = 1; % direct path
            channel_taps(rand_idxs) = channel_vals(:);
        
            y3(:,ii) = conv(y2(:,ii),channel_taps);
        end  
    else % no multipath, do nothing
        y3 = y2; 
    end
       
    %% add noise
    Nsamp = size(y3,1);
    y4 = y3 + sigman*(randn(Nsamp,num_rx) + 1j*randn(Nsamp,num_rx));
    
    %% time domain correlation for each unique rx pair
    for mm = 1:num_rx-1
        [corrs_t(nn,:,mm), lags] = xcorr(y3(:,mm+1),y3(:,1));
    end
end
% mean_delays_est = mean(delays_est);   % average delay estimate over trials
% std_delays_est = std(delays_est);     % standard deviation over trials
    
%% Plots
% Time domain - compare delayed version to non-delayed version
figure
subplot(2,1,1)
plot(x4+0.01,'b.-'); hold on
leg_str1{1} = 'Original';
for nn = 1:num_rx
    plot(real(y3(:,nn)), '.-')
    leg_str1{nn+1} = sprintf('%3.1f ns Delay',delays(nn));
end
xlabel('Sample Number')
ylabel('Amplitude')
title('Time Domain of Original and Delayed Signals')
legend(leg_str1{:})
axis tight

% plot the correlations
plot_idxs = 1:Ntrials;
plot_opts = {'b.-','r.-','g.-','c.-','y.-','m.-'};
subplot(212)
for nn = 1:length(plot_idxs)
    for mm = 1:num_rx-1
        plot(lags*Ts, abs(corrs_t(plot_idxs(nn),:,mm)), plot_opts{mm}, 'markersize', 10); hold on
        leg_str{mm} = sprintf('%3.1f ns Delay btwn Rx%i and Rx1',delays(mm+1),mm+1);
    end
end
ylims = ylim;
plot([0 0], [ylims(1) ylims(2)], 'k--')
title('Time Domain Correlations')
legend(leg_str{:},'location','best')
xlabel('Sample Lag')
ylabel('|XCorr|')
grid on
xb = 2*max(abs(delays));
axis([-xb xb -inf inf])
ticks = -xb:Ts/2:xb;
tickstep = 4;
ticklabels = cell(1,length(ticks));
ticks_c = num2cell(ticks);
ticklabels(1:tickstep:end) = ticks_c(1:tickstep:end);
 set(gca,'Xtick',ticks)
 set(gca,'XtickLabel',ticklabels)
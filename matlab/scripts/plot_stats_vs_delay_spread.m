% the number of paths includes the line-of-sight (LOS) path in this script
% for example if num_paths = 2, then there is the direct path and one
% non-LOS (NLOS) path
clear; close all
addpath('../functions')

plot_toa_countours = 0;
Ntrials = 100;            % keep at 1 for this plot, only for heat maps
% plot_hull = 1;          % plot convex hull or ref stations

%% Channel parameters
delay_spread = 300e-9;   % time difference between first received path and last (s)
multi_option = 3;   % determines the type of plot to generate
                    % 0 = no multipath
                    % 1 = two path with various delays between them
                    % 2 = varrying number of paths between 1 and max_num_paths
                    % with increasing delay spread
                    % 3 = two path with increasing delay spread
num_paths = 2;      % number of multipaths per reciever for option 1
multi_jump = 4;     % skip amount for option 1, 1:multi_jump:num_paths
num_delay_spreads = 10; % number of delay spreads to test in option 2
max_num_paths = inf; % max number paths for option 2

% Generate target emitter positions [x; y] (meters)
% numEmitters = 100;       % total number of random emitter locations to create
% targetPos = [randi([emitter_bounds(1) emitter_bounds(2)],1,numEmitters); ...    
%     randi([emitter_bounds(3) emitter_bounds(4)],1,numEmitters)]; 
% [Tx, Ty] = meshgrid(emitter_bounds(1):emitter_spacing:emitter_bounds(2), ...
%     emitter_bounds(3):emitter_spacing:emitter_bounds(4));

%% Receiver coords using rectangular coords 
% Equilateral Triangle
% a = 3.9624;     % length of one side of desired equilateral triangle
% a = 39.4908;
% b = sqrt(3)*a/2;
% refPos = [ 0  -a          -a/2; ...   % equilateral triangle
%            0  0      sqrt(3)*a/2];
% % center = [sum(refPos(1,:))/3; sum(refPos(2,:))/3];
% 
% center = [0;0];
% refPos = -center + refPos; % origin centered equilateral triangle
% % refPos = [[0;0] refPos];    % add a ref node at the center

%% Receiver coords using polar coords
a = 39.4908;
r1 = [22.8, -30*pi/180];
r2 = [22.8, -150*pi/180];
r3 = [22.8, 90*pi/180];
center = r1;
refPos = [ r1(1)*cos(r1(2)) r2(1)*cos(r2(2)) r3(1)*cos(r3(2)); ...   % equilateral triangle
           r1(1)*sin(r1(2)) r2(1)*sin(r2(2)) r3(1)*sin(r3(2))];
center = [center(1)*cos(center(2));center(1)*sin(center(2))];
refPos = -center + refPos; % origin centered equilateral triangle

%% Emitter coords rectangular coords
targetPos3 = [-3*a/4;sqrt(3)/4*a];   % chair
targetPos2 = [-a/2;0];               % base
targetPos4 = [-a/4;sqrt(3)/4*a];     % opp chair
targetPos1 = [sum(refPos(1,:))/3; sum(refPos(2,:))/3];     % center
targetPos5 = [-a/4+50;sqrt(3)/4*a+50];     % opp chair

% targetPos = [targetPos1 targetPos2 targetPos3 targetPos4 targetPos5];
% targetPos = [targetPos1 targetPos2 targetPos3 targetPos4];
targetPos = [targetPos1];
[numdims, numtargets] = size(targetPos);

%% Emitter coords polar coords
% targetPos3 = [-3*a/4;sqrt(3)/4*a];   % chair
% targetPos2 = [-a/2;0];               % base
% targetPos4 = [2*r1(1)^2-4*r1(1)*cos(120);120];     % opp chair
% targetPos1 = [sum(refPos(1,:))/3; sum(refPos(2,:))/3];     % center
% 
% targetPos = [targetPos1 targetPos2 targetPos3 targetPos4];
% % targetPos = [targetPos1];


%% The figures will be bounded by this region
bcenter = [sum(refPos(1,:))/3; sum(refPos(1,:))/3; sum(refPos(2,:))/3; sum(refPos(2,:))/3].';
% bounds = bcenter + [-3 3 -3 3;
%                   -3 3 -3 3;
%                   -3 3 -3 3
%                   -3 3 -3 3];

temp = ones(2,2*numtargets);
temp(1,:) = -temp(1,:);
temp2 = reshape(temp,4,[]).';
r1 = 30;
a = sqrt(2*r1(1)^2+4*r1(1)*cos(120*pi/180));
bounds = bcenter + 1.0*a*reshape(temp,4,[]).';
              
%% Emitter pulse properties
tx_pwr_dbm = -22;         % emitter transmit power in dBm (USRP max is 10 dBm)
fs_tx = 200e6/5;
Nsym = 10;              % number of symbols in signals
span = 10;              % total length of shaping filter in symbols
sps = 2;                % samples per symbol at the transmitter
fsym = fs_tx/sps;             % symbol rate of transmitter (signal bandwidth)
Tsym = sps/fs_tx;
beta = 0.4;             % excess bandwidth of tx pulse shaping filter
fc = 2.395e9;             % center frequency of transmitter

%% Receiver properties
fs = 200e6/2;                % receiver sample rates (Hz)
wlen = 2*ceil(fs/fsym)+1; % moving maximum window length in samples, odd number
nstds = 9;                % number of standard deviations to declare peak
percent_of_peak = 0.8;    % get the number of samples needed on either side 
                          % of correlation peaks for the peak value to drop 
                          % by this percent for use in super resolution
%% Needed by dpd
grid_bounds = [-43 3 -4 37];
grid_bounds = [-5 5 -5 5];
adder = 0;
grid_xmin = grid_bounds(1) + bcenter(1) - adder;
grid_xmax = grid_bounds(2) + bcenter(1) + adder;
grid_ymin = grid_bounds(3) + bcenter(3) - adder;
grid_ymax = grid_bounds(4) + bcenter(3) + adder;
grid_numx = 20;
grid_numy = 20;

grid_def = [grid_xmin grid_xmax;
            grid_ymin grid_ymax;
            grid_numx grid_numy];
        
%% General simulation properties
apply_calibration = 0; % for hardware sims
calib_path = [];
c = 299792458;          % speed of light m/s

tic
num_ds_samps = ceil(delay_spread*fs); % delay spread in rx samples count
if multi_option == 0
    multi_idxs = {0};
    num_jj = 1;
    num_kk = 1;
elseif multi_option == 1
    multi_idxs = num2cell(0:multi_jump:num_ds_samps);
    num_jj = 1;
    num_kk = length(multi_idxs);
elseif multi_option == 2
    delay_spread = linspace(0,delay_spread,num_delay_spreads);
    multi_idxs = {nan}; % this will be filled in randomly later
    num_jj = num_delay_spreads;
    num_kk = 1;
elseif multi_option == 3
    delay_spread = linspace(0,delay_spread,num_delay_spreads);
    multi_idxs = {nan}; % this will be filled in randomly later
    num_jj = num_delay_spreads;
    num_kk = 1;
else
    fprintf(1,'\n\nOption not implemented\n\n')
end
num_multi_idxs = length(multi_idxs);
initial_coords = [-1;1];
for ii = 1:numtargets
    for jj = 1:num_jj
        for kk = 1:num_kk
            [coords{ii,jj*kk}, bias_coords{ii,jj*kk}, covar_coords{ii,jj*kk}, ...
            mse_coords(ii,jj*kk), tdoas_true(:,ii), tdoas_coarse{ii,jj*kk}, ...
            tdoas_refined{ii,jj*kk}, prob_correlation(ii,jj*kk), ...
            prob_detection(ii,jj*kk), avg_snr_db{ii,jj*kk}, dpd_grid{ii,jj*kk}, ...
            multi_rands{ii,jj*kk}, unique(ii,jj*kk), corr_mag_sqs{ii,jj*kk}] = ...
            get_single_emitter2(targetPos(:,ii), refPos, Ntrials, tx_pwr_dbm, ...
            fc, fs, fsym, Nsym, span, sps, beta, wlen, nstds, percent_of_peak, ...
            apply_calibration, calib_path, grid_def, delay_spread(jj), ...
            num_paths, max_num_paths, multi_idxs{kk}, multi_option, ...
            initial_coords, 0);
        end
    end
end
toc

%% Plots
if multi_option == 1
temp = horzcat(bias_coords{:});
bx = temp(1,:);
by = temp(2,:);
temp = horzcat(covar_coords{:});
for jj = 1:num_multi_idxs
    temp = diag(covar_coords{1,jj});
    stdx(jj) = sqrt(temp(1));
    stdy(jj) = sqrt(temp(2));
end

% Plot multipath index vs statistics
figure
subplot(311)
plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), sqrt(mse_coords(1,:)), '.-'); hold all
ylims = ylim;
ylims(2) = min(60,ylims(2));
plot([Tsym:Tsym:delay_spread; Tsym:Tsym:delay_spread]/1e-9, [ylims(1) ylims(2)], 'k--');
title(sprintf('Two Paths with Increasing Delay Between Them, Tx Power %3.0f dBm', tx_pwr_dbm))
xlabel('Delay Between 1st and 2nd Path (ns)')
ylabel('RMSE (m)')
grid on
axis([-inf inf ylims(1) ylims(2)])

subplot(312)
plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), bx, '.-'); hold all
plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), by, '.-')
ylims = ylim;
ylims(1) = max(-60,ylims(1));
ylims(2) = min(60,ylims(2));
plot([Tsym:Tsym:delay_spread; Tsym:Tsym:delay_spread]/1e-9, [ylims(1) ylims(2)], 'k--');
xlabel('Delay Between 1st and 2nd Path (ns)')
ylabel('Bias (m)')
legend('b_x','b_y')
grid on
axis([-inf inf ylims(1) ylims(2)])

subplot(313)
plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), stdx, '.-'); hold all
plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), stdy, '.-')
ylims = ylim;
ylims(2) = min(60,ylims(2));
plot([Tsym:Tsym:delay_spread; Tsym:Tsym:delay_spread]/1e-9, [ylims(1) ylims(2)], 'k--');
xlabel('Delay Between 1st and 2nd Path (ns)')
ylabel('Standard Deviation (m)')
legend('\sigma_x','\sigma_y')
grid on
axis([-inf inf ylims(1) ylims(2)])

plot_opts = {'b.-', 'r.-', 'g.-', 'm.-'};
num_opts = length(plot_opts);
kk = 1;
figure
% title('Correlation Peaks across Path Delays')
for ii = 1:num_multi_idxs
	subplot(num_multi_idxs,1,ii)
    nsamps = length(corr_mag_sqs{ii});
    lags = (1:nsamps) - (nsamps+1)/2;
    plot(lags/(fs*1e-9), corr_mag_sqs{ii}, plot_opts{kk}); hold on
    axis([-3*delay_spread 3*delay_spread -inf inf]/1e-9)
    
    if kk == num_opts
        kk = 1;
    else
        kk = kk + 1;
    end
end
elseif multi_option == 2 || multi_option == 3
    
temp = horzcat(bias_coords{:});
bx = temp(1,:);
by = temp(2,:);
temp = horzcat(covar_coords{:});
for jj = 1:num_delay_spreads
    temp = diag(covar_coords{1,jj});
    stdx(jj) = sqrt(temp(1));
    stdy(jj) = sqrt(temp(2));
end

figure
subplot(311)
plot(delay_spread/1e-9, sqrt(mse_coords(1,:)), '.-'); hold all
ylims = ylim;
ylims(2) = min(60,ylims(2));
plot([Tsym:Tsym:max(delay_spread); Tsym:Tsym:max(delay_spread)]/1e-9, [ylims(1) ylims(2)], 'k--');
title(sprintf('Random Number of Paths and Path Delays, Tx Power %3.0f dBm', tx_pwr_dbm))
xlabel('Delay Spread (ns)')
ylabel('RMSE (m)')
grid on
axis([-inf inf ylims(1) ylims(2)])

subplot(312)
plot(delay_spread/1e-9, bx, '.-'); hold all
plot(delay_spread/1e-9, by, '.-')
ylims = ylim;
ylims(1) = max(-60,ylims(1));
ylims(2) = min(60,ylims(2));
plot([Tsym:Tsym:max(delay_spread); Tsym:Tsym:max(delay_spread)]/1e-9, [ylims(1) ylims(2)], 'k--');
xlabel('Delay Spread (ns)')
ylabel('Bias (m)')
legend('b_x','b_y')
grid on
axis([-inf inf ylims(1) ylims(2)])

subplot(313)
plot(delay_spread/1e-9, real(stdx), '.-'); hold all
plot(delay_spread/1e-9, real(stdy), '.-')
ylims = ylim;
ylims(2) = min(60,ylims(2));
plot([Tsym:Tsym:max(delay_spread); Tsym:Tsym:max(delay_spread)]/1e-9, [ylims(1) ylims(2)], 'k--');
xlabel('Delay Spread (ns)')
ylabel('Standard Deviation (m)')
legend('\sigma_x','\sigma_y')
grid on
axis([-inf inf ylims(1) ylims(2)])

plot_opts = {'b.-', 'r.-', 'g.-', 'm.-'};
num_opts = length(plot_opts);
kk = 1;
figure
for ii = 1:num_delay_spreads
	subplot(num_delay_spreads,1,ii)
    nsamps = length(corr_mag_sqs{ii});
    lags = (1:nsamps) - (nsamps+1)/2;
    plot(lags/(fs*1e-9), corr_mag_sqs{ii}, plot_opts{kk}); hold on
    axis([-1*delay_spread(end) 1*delay_spread(end) -inf inf]/1e-9)
    
    if kk == num_opts
        kk = 1;
    else
        kk = kk + 1;
    end
end
end


%% Plot localization results
figure
max_cols = 4;
num_rows = ceil(max(num_jj,num_kk)/max_cols);  % 4 plots across max
for kk = 1:max(num_jj,num_kk)
%     plot_idx = mod(kk,max_cols);
    subplot(num_rows,max_cols,kk)
    numrefs = size(refPos,2);
    for ii = 1:numrefs
        plot(refPos(1,ii), refPos(2,ii), 'ks', 'MarkerFaceColor', 'k', ...
            'MarkerSize',10, 'HandleVisibility','off'); 
        hold all
    end
    plot(refPos(1,1), refPos(2,1), 'ks', 'MarkerFaceColor', 'k')
    for ii = 1:numrefs
        h = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle',...
            'FontSize', 8);
        set(h, 'Color',[1, 1 ,1])
    end
    plot(targetPos(1),targetPos(2), 'kx', 'MarkerSize',18);

    % plot dpd grid
    plot(dpd_grid{1,1}{1}, dpd_grid{1,1}{2}, 'k.', 'MarkerFaceColor', 'k', ...
            'MarkerSize',6, 'HandleVisibility','off'); 
    
    for ii = 1:Ntrials
        plot(coords{1,kk}(1,ii),coords{1,kk}(2,ii), 'b.', 'MarkerSize',...
            12, 'HandleVisibility','off');
    end
    plot(coords{1,1}(1),coords{1,1}(2), 'b.', 'MarkerSize',12)
    plot(initial_coords(1), initial_coords(2), 'ko')
    axis equal
    axis(bounds(1,:))
    xlabel('x (m)')
    ylabel('y (m)')
    xlen = bounds(1,2) - bounds(1,1);
    ylen = bounds(1,4) - bounds(1,3);
    text(bounds(1,2)-0.35*xlen,bounds(1,4)-0.15*ylen, ...
    sprintf('b_x: %3.2f (m)\nb_y: %3.2f (m)\n\\sigma_x: %3.2e (m)\n\\sigma_y: %3.2e (m)', ...
    bias_coords{kk}(1), bias_coords{kk}(2), sqrt(covar_coords{kk}(1,1)), ...
    sqrt(covar_coords{kk}(2,2))), 'fontsize', 8)

    
    % Plot time contours
    if plot_toa_countours == 1
        tbounds = 2*bounds(kk,:);
        steps = .006*max(tbounds);
        [tx,ty] = meshgrid(tbounds(1):steps:tbounds(2),tbounds(3):steps:tbounds(4));
        tx2 = tx - targetPos(1,kk);
        ty2 = ty - targetPos(2,kk);
        tz = sqrt(tx2.^2+ty2.^2)/c*1e9;
        contour(h1,tx,ty,tz, 4:15, '--', 'showtext', 'on', 'HandleVisibility', 'off')
    end
end
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'Position', [0.5 0.9 0.1 0.1])

%% Plot tdoa estimation and statistics
figure
h2 = gca;
for kk = 1:num_kk
    xmin = inf;
    xmax = -inf;
    xmin = inf;
    xmax = -inf;
    ax(kk) = subplot(num_kk,1,kk);
    mpdelay = multi_idxs{kk}/fs/1e-9;
    for nn = 1:numrefs-1
        refined_tdoas = tdoas_refined{kk}(:,nn)*1e9;
        mean_tdoa = mean(refined_tdoas);
        std_tdoa = std(refined_tdoas);
        htrue = plot(tdoas_true, nn, 'kx', 'markersize', 8); hold all
        hcoarse = plot(tdoas_coarse{kk}(:,nn)*1e9, nn, 'ro');
        hrefine = plot(tdoas_refined{kk}(:,nn)*1e9, nn, 'b.');
        axis([-10 10 0.5 numrefs])
        xticks(-10:10)
%         xmin = min([xmin; tdoas_true; tdoas_refined{kk}(:,nn)*1e9]);
%         xmax = max([xmax; tdoas_true; tdoas_refined{kk}(:,nn)*1e9]);
%         axis([min(-1, 1.1*xmin) max(1.1*xmax,1) 0.5 numrefs])
        text(tdoas_true(nn), nn, sprintf('\\mu = %3.2f, \\sigma = %3.2f', mean_tdoa, std_tdoa), ...
            'horizontalalignment', 'center', 'verticalalignment', 'top',...
            'FontSize', 8);
    end
    title(sprintf('2nd Path Delay: %3.0f ns', mpdelay))
    y_values = 1:numrefs-1;
    ylabels = {'1,2', '1,3'};
    set(gca, 'Ytick', y_values, 'YTickLabel',ylabels);
    ylabel(sprintf('Receiver\n Correlation\n Pairs'))
end
xlabel('TDOA (ns)')
hleglines = [htrue(1) hcoarse(1) hrefine(1)];
legend(ax(1), hleglines, 'True TDOA', 'Coarse TDOA Estimate', 'Refined TDOA Estimate', 'Location', 'Northeast')
set(gcf, 'Position',  [100, 100, 1300, 800])

% save('tdoa_data.mat', 'bx','by','stdx','stdy','coords','delay_spread', ...
%     'bias_coords','covar_coords','refPos','targetPos','num_jj','num_kk',...
%     'multi_option', 'mse_coords','num_delay_spreads', 'Tsym','tx_pwr_dbm',...
%     'corr_mag_sqs','fs','Ntrials','initial_coords','bounds','avg_snr_db')

save('dpd_data.mat', 'bx','by','stdx','stdy','coords','delay_spread', ...
    'bias_coords','covar_coords','refPos','targetPos','num_jj','num_kk',...
    'multi_option', 'mse_coords','num_delay_spreads', 'Tsym','tx_pwr_dbm',...
    'corr_mag_sqs','fs','Ntrials','initial_coords','bounds','avg_snr_db','grid_numx',...
    'grid_numy','dpd_grid')
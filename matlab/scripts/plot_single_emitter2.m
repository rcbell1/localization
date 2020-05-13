clear; close all
addpath('../functions')

show_plots = 0;         % show plots for debugging
show_circles = 0;       % plot circles centered on emitter to visualize tdoa
show_hyperbolas = 0;    % plot hyperbolas to visualize intersection point
Ntrials = 10;            % keep at 1 for this plot, only for heat maps

%% Channel parameters
delay_spread = 500e-9;   % time difference between first received path and last (s)
multi_option = 0;   % determines the type of plot to generate
                    % 0 = no multipath
                    % 1 = two path with various delays between them
                    % 2 = varrying number of paths between 1 and num_paths
                    % with increasing delay spread
num_paths = 2;      % number of multipaths per reciever for option 1
multi_jump = 5;     % skip amount for option 1, 1:multi_jump:num_paths
num_delay_spreads = 10; % number of delay spreads to test in option 2
max_num_paths = inf; % max number paths for option 2

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

%% Emiiter coords
targetPos = [sum(refPos(1,:))/3; sum(refPos(2,:))/3];     % center

%% The figures will be bounded by this region
bcenter = [sum(refPos(1,:))/3; sum(refPos(1,:))/3; sum(refPos(2,:))/3; sum(refPos(2,:))/3].';
% bounds = bcenter + [-3 3 -3 3;
%                   -3 3 -3 3;
%                   -3 3 -3 3
%                   -3 3 -3 3];
numtargets = 1;
temp = ones(2,2*numtargets);
temp(1,:) = -temp(1,:);
temp2 = reshape(temp,4,[]).';
r1 = 30;
a = sqrt(2*r1(1)^2+4*r1(1)*cos(120*pi/180));
bounds = bcenter + 1.0*a*reshape(temp,4,[]).';

%% Emitter pulse properties
tx_pwr_dbm = -20;         % emitter transmit power in dBm (USRP max is 10 dBm)
fs_tx = 200e6/10;
Nsym = 1000;              % number of symbols in signals
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

%% Estimate the location
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
else
    fprintf(1,'\n\nOption not implemented\n\n')
end
num_multi_idxs = length(multi_idxs);
initial_coords = [-41;1];

% [coords, bias_coords, covar_coords, mse_coords, tdoas] = ...
%     get_single_emitter2(targetPos, refPos, 1, tx_pwr_dbm, fc, ...
%     fs, fsym, Nsym, span, sps, beta, wlen, nstds, percent_of_peak,...
%     show_plots);

[coords, bias_coords, covar_coords, mse_coords, tdoas_true, tdoas_coarse, ...
tdoas_refined, prob_correlation, prob_detection, dpd_grid, a, unique, ...
corr_mag_sq] = get_single_emitter2(targetPos, refPos, Ntrials, tx_pwr_dbm, ...
fc, fs, fsym, Nsym, span, sps, beta, wlen, nstds, percent_of_peak, ...
apply_calibration, calib_path, grid_def, delay_spread, num_paths, ...
max_num_paths, multi_idxs, multi_option, initial_coords, 0);

toc

%% Plots
% Plot localization results
figure
subplot(1,2,1)
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
plot(targetPos(1),targetPos(2), 'bo', 'MarkerSize',8);
for ii = 1:size(coords,2)
    plot(coords(1,ii),coords(2,ii), 'b.', 'MarkerSize',12, 'HandleVisibility','off');
end
plot(coords(1,1),coords(2,1), 'b.', 'MarkerSize',12)
axis equal
axis(bounds)
xlabel('x (m)')
ylabel('y (m)')
title('Localization Results')
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'Position', [0.25 0.06 0.1 0.1])
% grid on

% Plot time contours
tbounds = 2*bounds;
steps = round(.1*max(tbounds));
[tx,ty] = meshgrid(tbounds(1):steps:tbounds(2),tbounds(3):steps:tbounds(4));
tx2 = tx - targetPos(1);
ty2 = ty - targetPos(2);
tz = sqrt(tx2.^2+ty2.^2)/c*1e9;
contour(tx,ty,tz, '--', 'showtext', 'on', 'HandleVisibility', 'off')

subplot(1,2,2)
for ii = 1:numrefs
    plot(refPos(1,ii), refPos(2,ii), 'rx','HandleVisibility','off', 'MarkerSize',14); 
    hold all
    text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
end
plot(refPos(1,1), refPos(2,1), 'rx')
plot(targetPos(1),targetPos(2), 'bo', 'MarkerSize',8);
for ii = 1:size(coords,2)
    plot(coords(1,ii),coords(2,ii), 'b.', 'MarkerSize',12, 'HandleVisibility','off');
end
plot(coords(1,1),coords(2,1), 'b.', 'MarkerSize',12)
axis equal
axis(bounds)
xlabel('x (m)')
ylabel('y (m)')
title('Localization Results')
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'Position', [0.69 0.06 0.1 0.1])
grid on

if show_circles == 1
    viscircles(repmat(targetPos, 1, numrefs)',vecnorm(targetPos-refPos,2,1),...
    'color', 'b', 'linestyle', '--', 'linewidth', 0.5);
end

if show_hyperbolas == 1
    plot_multiple_hyperbolas(refPos, targetPos, tdoas*c, bounds)
end

figure
xmin = inf;
xmax = -inf;
for nn = 1:numrefs-1
    refined_tdoas = tdoas_refined(:,nn)*1e9;
    mean_tdoa = mean(refined_tdoas);
    std_tdoa = std(refined_tdoas);
    htrue = plot(tdoas_true(nn), nn, 'kx', 'markersize', 8); hold all
    hcoarse = plot(tdoas_coarse(:,nn)*1e9, nn, 'ro');
    hrefine = plot(tdoas_refined(:,nn)*1e9, nn, 'b.');
    xmin = min([xmin; tdoas_true(nn); tdoas_refined(:,nn)*1e9]);
    xmax = max([xmax; tdoas_true(nn); tdoas_refined(:,nn)*1e9]);
    axis([min(-1, 1.1*xmin) max(1.1*xmax,1) 0.5 numrefs])
    text(mean_tdoa, nn, sprintf('\\mu = %3.2f, \\sigma = %3.2f', mean_tdoa, std_tdoa), ...
        'horizontalalignment', 'center', 'verticalalignment', 'top',...
        'FontSize', 8);
end
y_values = 1:numrefs-1;
ylabels = {'1,2', '1,3'};
set(gca, 'Ytick', y_values, 'YTickLabel',ylabels);
xlabel('TDOA (ns)')
ylabel('Receiver Correlation Pairs')
hleglines = [htrue(1) hcoarse(1) hrefine(1)];
legend(hleglines, 'True TDOA', 'Coarse TDOA Estimate', 'Refined TDOA Estimate', 'Location', 'North')

set(gcf, 'Position',  [100, 100, 1300, 800])





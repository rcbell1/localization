clear; close all
addpath('../functions')

show_plots = 1;         % show plots for debugging
show_circles = 1;       % plot circles centered on emitter to visualize tdoa
show_hyperbolas = 1;    % plot hyperbolas to visualize intersection point
% targetPos = [0; -5];    % target position (meters)
targetPos = [10000; 10000];
% refPos = [-5 5 0 ; ...  % reference receiver positions [x; y] (meters)
%           -5 -5 5 ]; 
% refPos = [-50 50 0; ... % triangle
%           -50 -50 50];
% refPos = [50 -50 -50 50 ; ... % square
%           -50 -50 50 50];
a = 5;
refPos = -[a/2; a/(3*sqrt(2))] + ... % centered triangle
    [ 0  a    a/2; ... 
      0  0      a/sqrt(2)];
refPos = [[0;0] refPos];
% refPos = [0  -40  40 -70 70; ... % 5-pnt star
%           0 -40 -40  40 40];
% refPos = [0  -40  40 -70 70 100 -100;
%            0 -40 -40  40 40 -100 100];
% refPos = refPos - refPos(:,1);
% targetPos = [60; 60] - refPos(:,1);
bounds = [-10000 10000 -10000 10000];
% bounds = [-150 150 -150 150];

% Emitter pulse properties
tx_pwr_dbm = 7;         % emitter transmit power in dBm
fc = 915e6;             % center frequency of transmitter
span = 10;              % total length of shaping filter in symbols
sps = 4;                % samples per symbol at the receiver sample rate
beta = 0.4;             % excess bandwidth of tx pulse shaping filter
Nsym = 40;              % number of symbols in signals
fsym = 9e6;             % symbol rate of transmitter (signal bandwidth)

% Receiver properties
fs = 20e6;                % receiver sample rates (Hz)
wlen = 2*ceil(fs/fsym)+1; % moving maximum window length in samples
nstds = 12;                % number of standard deviations to declare peak
percent_of_peak = 0.8;    % get the number of samples needed on either side 
                          % of correlation peaks for the peak value to drop 
                          % by this percent for use in super resolution

% General simulation properties
c = 299792458;          % speed of light m/s

%% Estimate the location
[coords, bias_coords, covar_coords, mse_coords, tdoas] = ...
    get_single_emitter2(targetPos, refPos, 1, tx_pwr_dbm, fc, ...
    fs, fsym, Nsym, span, sps, beta, wlen, nstds, percent_of_peak,...
    show_plots);

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
steps = round(.1*tbounds);
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







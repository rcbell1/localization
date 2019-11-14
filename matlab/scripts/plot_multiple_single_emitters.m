clear; close all
addpath('../functions')

plot_hull = 1;          % plot convex hull or ref stations
emitter_bounds = [-40 40 -40 70];   % bounds of emitter locations
plot_bounds = [-250 250 -250 250];  % outer boundary of plot figure

% Generate target emitter positions [x; y] (meters)
numEmitters = 4;       % total number of random emitter locations to create
targetPos = [randi([emitter_bounds(1) emitter_bounds(2)],1,numEmitters); ...    
    randi([emitter_bounds(3) emitter_bounds(4)],1,numEmitters)]; 
[Tx, Ty] = meshgrid(emitter_bounds(1):emitter_spacing:emitter_bounds(2), ...
    emitter_bounds(3):emitter_spacing:emitter_bounds(4));

% Generate reference station positions [x; y] (meters)
% refPos = [-40 40 0 -70 70 -20 20 -35 -8 22; ... 
%           -40 -40 80 40 40 40 40 16 -7 11]; 
refPos = [-50 50 0; ... % triangle
          -50 -50 50];

tx_pwr_dbm = -10;       % emitter transmit power in dBm
fc = 915e6;             % center frequency of transmitter
fs = 20e6;              % receiver sample rates (Hz)
Nsym = 40;              % number of symbols in signals

span = 20;              % total length of shaping filter in symbols
sps = 5;                % samples per symbol at the receiver sample rate
fhigh = 2500*fs;           % high speed sample rate where delays are added (Hz)
wlen = 2*sps+1;           % moving maximum window length in samples
nstds = 6;              % number of standard deviations to declare peak

Ntrials = 1;            % keep at 1 for this plot, only for heat maps
show_plots = 0;         % don't use this for this script, too many plots

tic
for ii = 1:size(targetPos,2)
    [coords(:,ii), bias_coords(:,ii), ~, mse_coords(ii)] = ...
        get_single_emitter(targetPos(:,ii), refPos, Ntrials, ...
        tx_pwr_dbm, fc, fs, Nsym, span, sps, fhigh, wlen, nstds, show_plots);
end
toc

figure
numrefs = size(refPos,2);
for ii = 1:numrefs
    plot(refPos(1,ii), refPos(2,ii), 'r*','HandleVisibility','off', 'MarkerSize',8); 
    hold all
    text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
end
plot(refPos(1,1), refPos(2,1), 'r*')

for ii = 1:size(coords,2)
    plot(targetPos(1,ii),targetPos(2,ii), 'bo', 'MarkerSize',8, 'HandleVisibility','off');
    plot(coords(1,ii),coords(2,ii), 'b.', 'MarkerSize',12, 'HandleVisibility','off');
end
plot(targetPos(1,1),targetPos(2,1), 'bo')
plot(coords(1,1),coords(2,1), 'b.', 'MarkerSize',12)
axis equal
axis(plot_bounds)
xlabel('x (m)')
ylabel('y (m)')
title('Localization Results')
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'location', 'best')
grid on

if plot_hull == 1
    k = convhull(refPos(1,:),refPos(2,:));
end
plot(refPos(1,k),refPos(2,k),'r', 'HandleVisibility','off')

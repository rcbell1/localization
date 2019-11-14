clear; close all
addpath('../functions')
addpath('../../intersecting_hyperbolas')

show_plots = 1;         % show plots for debugging
show_circles = 1;       % plot circles centered on emitter to visualize tdoa
show_hyperbolas = 1;    % plot hyperbolas to visualize intersection point
targetPos = [0; 0];     % target position (meters)
refPos = [-50 50 0 ; ... % reference receiver positions [x; y] (meters)
          -50 -50 50 ];  
% refPos = [0  -40  40 -70 70; ... % 5-pnt star
%           80 -40 -40  40 40];
tx_pwr_dbm = -90;       % emitter transmit power in dBm
Ntrials = 1;            % number of simulations per emitter location
fs = 20e6;              % receiver sample rates (Hz)
Nsym = 40;              % number of symbols in signals

fc = 915e6;             % center frequency of transmitter
span = 20;              % total length of shaping filter in symbols
sps = 5;                % samples per symbol at the receiver sample rate
fhigh = 2500*fs;        % high speed sample rate where delays are added (Hz)
wlen = 2*sps+1;         % moving maximum window length in samples
nstds = 6;              % number of standard deviations to declare peak
c = 299792458;          % speed of light m/s

bounds = [-150 150 -150 150];

[coords, bias_coords, covar_coords, mse_coords, tdoas] = ...
    get_single_emitter(targetPos, refPos, Ntrials, tx_pwr_dbm, fc, ...
    fs, Nsym, span, sps, fhigh, wlen, nstds, show_plots);

figure
subplot(1,2,1)
numrefs = size(refPos,2);
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
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'location', 'best')
grid on

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
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'location', 'best')
grid on

if show_circles == 1
    viscircles(repmat(targetPos, 1, numrefs)',vecnorm(targetPos-refPos,2,1),...
    'color', 'b', 'linestyle', '--', 'linewidth', 0.5);
end

if show_hyperbolas == 1
    plot_multiple_hyperbolas(refPos, targetPos, tdoas*c, bounds)
end


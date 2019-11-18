clear; close all
addpath('../functions')

show_plots = 0;         % show plots for debugging
show_circles = 1;       % plot circles centered on emitter to visualize tdoa
show_hyperbolas = 1;    % plot hyperbolas to visualize intersection point
targetPos = [-7; 0];     % target position (meters)
refPos = [-5 5 0 ; ... % reference receiver positions [x; y] (meters)
          -5 -5 5 ];  
% refPos = [0  -40  40 -70 70; ... % 5-pnt star
%           80 -40 -40  40 40];
bounds = [-15 15 -15 15];
tx_pwr_dbm = 100;       % emitter transmit power in dBm
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

[coords, bias_coords, covar_coords, mse_coords, tdoas] = ...
    get_single_emitter(targetPos, refPos, Ntrials, tx_pwr_dbm, fc, ...
    fs, Nsym, span, sps, fhigh, wlen, nstds, show_plots);

%% Plots
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
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'location', 'best')
% grid on

% Plot time contours
tbounds = 2*bounds;
[tx,ty] = meshgrid(tbounds(1):tbounds(2),tbounds(3):tbounds(4));
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
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'location', 'best')
grid on

if show_circles == 1
    viscircles(repmat(targetPos, 1, numrefs)',vecnorm(targetPos-refPos,2,1),...
    'color', 'b', 'linestyle', '--', 'linewidth', 0.5);
end

if show_hyperbolas == 1
    plot_multiple_hyperbolas(refPos, targetPos, tdoas*c, bounds)
end






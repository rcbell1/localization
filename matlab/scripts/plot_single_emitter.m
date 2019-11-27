clear; close all
addpath('../functions')

show_plots = 1;         % show plots for debugging
show_circles = 1;       % plot circles centered on emitter to visualize tdoa
show_hyperbolas = 1;    % plot hyperbolas to visualize intersection point
targetPos = [-5; 25];    % target position (meters)
% refPos = [-5 5 0 ; ...  % reference receiver positions [x; y] (meters)
%           -5 -5 5 ]; 
% refPos = [-50 50 0; ... % triangle
%           -50 -50 50];
% refPos = [50 -50 -50 50 ; ... % square
%           -50 -50 50 50];
a = 50;
refPos = -[a/2; a/(3*sqrt(2))] + ... % centered triangle
    [ 0  a    a/2; ... 
      0  0      a/sqrt(2)];
refPos = [[0;0] refPos]
% refPos = [0  -40  40 -70 70; ... % 5-pnt star
%           0 -40 -40  40 40];
% refPos = [0  -40  40 -70 70 100 -100;
%            0 -40 -40  40 40 -100 100];
% refPos = refPos - refPos(:,1);
% targetPos = [60; 60] - refPos(:,1);
% bounds = [-15 15 -15 15];
bounds = [-150 150 -150 150];

% Emitter pulse properties
tx_pwr_dbm = 0;       % emitter transmit power in dBm
fc = 915e6;             % center frequency of transmitter
span = 10;              % total length of shaping filter in symbols
sps = 4;                % samples per symbol at the receiver sample rate
beta = 0.4;             % excess bandwidth of tx pulse shaping filter
Nsym = 40;              % number of symbols in signals
fsym = 1e6;             % symbol rate of transmitter (signal bandwidth)

% Receiver properties
fs = 50e6;                % receiver sample rates (Hz)
wlen = 2*ceil(fs/fsym)+1; % moving maximum window length in samples
nstds = 3;                % number of standard deviations to declare peak

% General simulation properties
Ntrials = 1;            % number of simulations per emitter location
c = 299792458;          % speed of light m/s
fhigh = 50e9;           % high speed sample rate where delays are added (Hz)

% Determine the value of fhigh that is above requested and an integer
% multiple of fs and fsym
ftemp = lcm(fs,fsym);
if ftemp < fhigh
    mult = ceil(fhigh/ftemp);
    fhigh = ftemp*mult; % nearest integer multiple of fs and fsym above fhigh
end

[coords, bias_coords, covar_coords, mse_coords, tdoas] = ...
    get_single_emitter(targetPos, refPos, Ntrials, tx_pwr_dbm, fc, ...
    fs, fsym, Nsym, span, sps, beta, fhigh, wlen, nstds, show_plots);

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

% subplot(2,2,3)
% for ii = 2:numrefs
%     pairStrings(ii-1) = sprintf("%i1",ii);
% end
% uitable('Station Pair',pairStrings)






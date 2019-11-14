clear; close all
addpath('../functions')

plot_hull = 1;          % plot convex hull or ref stations
emitter_bounds = [-150 150 -150 150];   % bounds of emitter locations
% emitter_bounds = [-700 700 -700 700];   % bounds of emitter locations
emitter_spacing = 2;   % spacing between test locations in meters
tx_pwr_dbm = -50;       % emitter transmit power in dBm
Ntrials = 100;          % number of simulations per emitter location
fc = 915e6;             % center frequency of transmitter
[Tx, Ty] = meshgrid(emitter_bounds(1):emitter_spacing:emitter_bounds(2), ...
    emitter_bounds(3):emitter_spacing:emitter_bounds(4));
% Reference receiver positions [x; y] (meters)
% refPos = [-40 40 0 -70 70 -20 20 -35 -8 22; ... % 10-pnt star
%           -40 -40 80 40 40 40 40 16 -7 11];   
refPos = [-50 50 0; ... % triangle
          -50 -50 50]; 
% refPos = [50 -50 -50 50 ; ... % square
%           -50 -50 50 50]; 
% refPos = [0  -40  40 -70 70; ... % 5-pnt star
%           80 -40 -40  40 40]; 
fs = 20e6;              % receiver sample rates (Hz)
Nsym = 40;              % number of symbols in signals

span = 20;              % total length of shaping filter in symbols
sps = 5;                % samples per symbol at the receiver sample rate
fhigh = 10e9;           % high speed sample rate where delays are added (Hz)
wlen = 2*sps+1;         % moving maximum window length in samples, odd number
nstds = 3;              % number of standard deviations to declare peak

show_plots = 0;         % leave this off for this sim, too many plots
tic
[nrows, ncols] = size(Tx);
parfor ii = 1:nrows
    for jj = 1:ncols
        targetPos = [Tx(ii,jj); Ty(ii,jj)];
        [~, ~, ~, mse_coords(ii,jj), ~] = ...
            get_single_emitter(targetPos, refPos, Ntrials, tx_pwr_dbm, ...
            fc, fs, Nsym, span, sps, fhigh, wlen, nstds, show_plots);
    end
end
run_time = toc/60;

figure
imagesc([emitter_bounds(1) emitter_bounds(2)],[emitter_bounds(3) ...
    emitter_bounds(4)], sqrt(mse_coords)); hold all
h = colorbar;           % displays the colorbar legend
colormap('jet')         % sets the colors used
caxis([0 40])           % sets the limits of the colormap
set(gca,'YDir','normal')% keeps y-axis correct orientation
xlabel('x (m)')
ylabel('y (m)')
title('Localization RMSE Performance')
title(h,sprintf('RMSE\n(m)'))

numrefs = size(refPos,2);
for ii = 1:numrefs
    plot(refPos(1,ii), refPos(2,ii), 'k.','HandleVisibility','off', ...
        'MarkerFaceColor', [0 1 0], 'MarkerSize',28);
    h = text(refPos(1,ii), refPos(2,ii)+3, sprintf('%i', ii), ...
        'horizontalalignment', 'center', 'verticalalignment', 'top');
    set(h, 'Color',[1, 1 ,1])
end

h = text(emitter_bounds(1),0.9*emitter_bounds(4), ...
    sprintf('Tx EIRP: %2.0f dBm', tx_pwr_dbm));
set(h, 'Color',[1, 1 ,1])

% Check condition number of A to make sure its not a problem. This only
% depends on the reference station placements
x1 = refPos(1,1);
y1 = refPos(2,1);
for ii = 2:size(refPos,2)
    A(:,ii-1) = refPos(:,1) - refPos(:,ii);    
end
A = A.';

fprintf(1,'\nTotal Runtime: %2.1f min \nMatrix Condition Number: %4.2f\n\n', ...
    run_time, cond(A));
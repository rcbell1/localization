clear; close all
addpath('../functions')

%% General simulation properties
% emitter_bounds = [-150 150 -150 150];   % bounds of emitter locations (m)
emitter_bounds = [-100 100 -100 100];   % bounds of emitter locations
emitter_spacing = 2;   % spacing between test emitter locations in meters

% Reference receiver positions [x; y] (meters)
% refPos = [-40 40 0 -70 70 -20 20 -35 -8 22; ... % 10-pnt star
%           -40 -40 80 40 40 40 40 16 -7 11];   
% refPos = [-50 50 0; ... % triangle
%           -50 -50 50]; 
% refPos = [50 -50 -50 50 ; ... % square
%           -50 -50 50 50]; 
a = 50;
refPos = -[a/2; a/(3*sqrt(2))] + ... % centered equilateral triangle
    [ 0  a    a/2; ... 
      0  0      a/sqrt(2)];
refPos = [[0;0] refPos];
% refPos = [0 -50 50 0; ... % centered triangle
%           0 -50 -50 50];
% refPos = [0  -40  40 -70 70; ... % 5-pnt star
%           80 -40 -40  40 40]; 
Ntrials = 100;          % number of noise instances per emitter location
fhigh = 50e9;           % high speed sample rate where delays are added (Hz)

%% Emitter pulse properties
tx_pwr_dbm = 0;         % emitter transmit power in dBm
Nsym = 40;              % number of symbols in signals
fsym = 4e6;             % symbol rate of transmitter (signal bandwidth)
span = 10;              % total length of shaping filter in symbols
sps = 4;                % samples per symbol at the transmitter
beta = 0.4;             % excess bandwidth of tx pulse shaping filter
fc = 50e9;             % center frequency of transmitter

%% Receiver properties
fs = 20e6;                % receiver sample rates (Hz)
wlen = 2*ceil(fs/fsym)+1; % moving maximum window length in samples, odd number
nstds = 3;                % number of standard deviations to declare peak

%% Simulation
% Determine the value of fhigh that is above requested and an integer
% multiple of fs and fsym
ftemp = lcm(fs,fsym);
if ftemp < fhigh
    mult = ceil(fhigh/ftemp);
    fhigh = ftemp*mult; % nearest integer multiple of fs and fsym above fhigh
end

% Create the grid of points based on the emitter bounds set
[Tx, Ty] = meshgrid(emitter_bounds(1):emitter_spacing:emitter_bounds(2), ...
    emitter_bounds(3):emitter_spacing:emitter_bounds(4));

tic
[nrows, ncols] = size(Tx);
parfor ii = 1:nrows
    for jj = 1:ncols
        targetPos = [Tx(ii,jj); Ty(ii,jj)];
        [~, ~, ~, mse_coords(ii,jj), ~, unique(ii,jj)] = ...
            get_single_emitter(targetPos, refPos, Ntrials, tx_pwr_dbm, ...
            fc, fs, fsym, Nsym, span, sps, beta, fhigh, wlen, nstds, ...
            0);
    end
end
run_time = toc/60;

%% Plots
% First plot the RMSE
figure
imagesc([emitter_bounds(1) emitter_bounds(2)],[emitter_bounds(3) ...
    emitter_bounds(4)], sqrt(mse_coords)); hold all
h = colorbar;           % displays the colorbar legend
colormap('jet')         % sets the colors used
caxis([0 2])           % sets the limits of the colormap
set(gca,'YDir','normal')% keeps y-axis correct orientation
xlabel('x (m)')
ylabel('y (m)')
title('Localization RMSE Performance')
title(h,sprintf('RMSE\n(m)'))

numrefs = size(refPos,2);
for ii = 1:numrefs
    plot(refPos(1,ii), refPos(2,ii), 'k.','HandleVisibility','off', ...
        'MarkerFaceColor', [0 1 0], 'MarkerSize',28);
    h = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
        'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
        'FontSize', 6);
    set(h, 'Color',[1, 1 ,1])
end

h = text(emitter_bounds(1),0.9*emitter_bounds(4), ...
    sprintf('Tx EIRP: %2.0f dBm', tx_pwr_dbm));
set(h, 'Color',[1, 1 ,1])

% Plot the uniqueness of solutions 
% figure
% imagesc([emitter_bounds(1) emitter_bounds(2)],[emitter_bounds(3) ...
%     emitter_bounds(4)], unique); hold all
% h = colorbar;           % displays the colorbar legend
% colormap('autumn')         % sets the colors used
% % caxis([0 40])           % sets the limits of the colormap
% set(gca,'YDir','normal')% keeps y-axis correct orientation
% xlabel('x (m)')
% ylabel('y (m)')
% title('Solution Uniqueness')
% % title(h,sprintf('RMSE\n(m)'))
% 
% numrefs = size(refPos,2);
% for ii = 1:numrefs
%     plot(refPos(1,ii), refPos(2,ii), 'k.','HandleVisibility','off', ...
%         'MarkerFaceColor', [0 1 0], 'MarkerSize',28);
%     h = text(refPos(1,ii), refPos(2,ii)+3, sprintf('%i', ii), ...
%         'horizontalalignment', 'center', 'verticalalignment', 'top');
%     set(h, 'Color',[1, 1 ,1])
% end


% Check condition number of A to make sure its not a problem. This only
% depends on the reference station placements
% x1 = refPos(1,1);
% y1 = refPos(2,1);
% for ii = 2:size(refPos,2)
%     A(:,ii-1) = refPos(:,1) - refPos(:,ii);    
% end
% A = A.';
% fprintf(1,'\nTotal Runtime: %2.1f min \nMatrix Condition Number: %4.2f\n\n', ...
%     run_time, cond(A));

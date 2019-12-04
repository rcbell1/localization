clear; close all
addpath('../functions')

%% General simulation properties
Ntrials = 100;          % number of noise instances per emitter location
% emitter_bounds = [-10000 10000 -10000 10000];   % Good for prob detection plots
emitter_bounds = [-10 10 -10 10];   % bounds of emitter locations
emitter_spacing = 0.1;   % spacing between test emitter locations in meters

% Reference receiver positions [x; y] (meters)
% refPos = [-40 40 0 -70 70 -20 20 -35 -8 22; ... % 10-pnt star
%           -40 -40 80 40 40 40 40 16 -7 11];   
% refPos = [-50 50 0; ... % triangle
%           -50 -50 50]; 
% refPos = [50 -50 -50 50 ; ... % square
%           -50 -50 50 50]; 
a = 5;     % length of one side of desired equilateral triangle
b = sqrt(3)*a/2;
refPos = [ 0  a          a/2; ...   % equilateral triangle
           0  0      sqrt(3)*a/2];
center = [sum(refPos(1,:))/3; sum(refPos(2,:))/3];
refPos = -center + refPos; % origin centered equilateral triangle
refPos = [[0;0] refPos];
% refPos = [0 -50 50 0; ... % centered triangle
%           0 -50 -50 50];
% refPos = [0  -40  40 -70 70; ... % 5-pnt star
%           80 -40 -40  40 40]; 

%% Emitter pulse properties
tx_pwr_dbm = 7;         % emitter transmit power in dBm (USRP max is 10 dBm)
Nsym = 40;              % number of symbols in signals
fsym = 5e6;             % symbol rate of transmitter (signal bandwidth)
span = 10;              % total length of shaping filter in symbols
sps = 4;                % samples per symbol at the transmitter
beta = 0.4;             % excess bandwidth of tx pulse shaping filter
fc = 915e6;             % center frequency of transmitter

%% Receiver properties
fs = 20e6;                % receiver sample rates (Hz)
wlen = 2*ceil(fs/fsym)+1; % moving maximum window length in samples, odd number
nstds = 9;                % number of standard deviations to declare peak
percent_of_peak = 0.8;    % get the number of samples needed on either side 
                          % of correlation peaks for the peak value to drop 
                          % by this percent for use in super resolution

%% Simulation
% Create the grid of points based on the emitter bounds set
[Tx, Ty] = meshgrid(emitter_bounds(1):emitter_spacing:emitter_bounds(2), ...
    emitter_bounds(3):emitter_spacing:emitter_bounds(4));

tic
[nrows, ncols] = size(Tx);
parfor ii = 1:nrows
    for jj = 1:ncols
        targetPos = [Tx(ii,jj); Ty(ii,jj)];
        [~, ~, ~, mse_coords(ii,jj), ~, prob_correlation(ii,jj), ...
            prob_detection(ii,jj), unique(ii,jj)] = ...
            get_single_emitter2(targetPos, refPos, Ntrials, tx_pwr_dbm, ...
            fc, fs, fsym, Nsym, span, sps, beta, wlen, nstds, ...
            percent_of_peak, 0);
    end
end
run_time = toc/60;

%% Plots
rmse_coords = sqrt(mse_coords);
[ndims, numrefs] = size(refPos);

% First plot the RMSE using pcolor
figure
hpc = pcolor(Tx, Ty, rmse_coords); hold all
colormap('jet')
h = colorbar;
caxis([0 4])           % sets the limits of the colormap
set(gca,'YDir','normal','color', 'w')% keeps y-axis correct orientation
set(hpc, 'EdgeColor', 'none')
hpc.FaceColor = 'interp';
xlabel('x (m)')
ylabel('y (m)')
title('Localization RMSE Performance')
title(h, sprintf('RMSE\n(m)'))

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

% Plot probability of detection using pcolor
figure
hpc = pcolor(Tx, Ty, prob_detection); hold all
colormap('jet')
h = colorbar;
caxis([0 1])           % sets the limits of the colormap
set(gca,'YDir','normal','color', 'w')% keeps y-axis correct orientation
set(hpc, 'EdgeColor', 'none')
xlabel('x (m)')
ylabel('y (m)')
title('Probability of Detection')
title(h,sprintf('Probability'))

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

% Plot probability of correlation using pcolor
figure
hpc = pcolor(Tx, Ty, prob_correlation); hold all
colormap('jet')
h = colorbar;
caxis([0 1])           % sets the limits of the colormap
set(gca,'YDir','normal','color', 'w')% keeps y-axis correct orientation
set(hpc, 'EdgeColor', 'none')
xlabel('x (m)')
ylabel('y (m)')
title('Probability of Correlation')
title(h,sprintf('Probability'))

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










% Plot the RMSE using imagesc
% I can't find a way to handle nans that doesn't effect other colors using
% imagesc. This isn't a problem when using pcolor
% figure
% cmap = colormap('jet');         % sets the colors used
% colormap(cmap)
% imgHand = imagesc([emitter_bounds(1) emitter_bounds(2)],[emitter_bounds(3) ...
%     emitter_bounds(4)], rmse_coords); hold all
% h = colorbar;           % displays the colorbar legend
% 
% caxis([0 4])           % sets the limits of the colormap
% set(gca,'YDir','normal')% keeps y-axis correct orientation
% xlabel('x (m)')
% ylabel('y (m)')
% title('Localization RMSE Performance')
% title(h,sprintf('RMSE\n(m)'))
% 
% for ii = 1:numrefs
%     plot(refPos(1,ii), refPos(2,ii), 'k.','HandleVisibility','off', ...
%         'MarkerFaceColor', [0 1 0], 'MarkerSize',28);
%     h = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
%         'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
%         'FontSize', 6);
%     set(h, 'Color',[1, 1 ,1])
% end
% h = text(emitter_bounds(1),0.9*emitter_bounds(4), ...
%     sprintf('Tx EIRP: %2.0f dBm', tx_pwr_dbm));
% set(h, 'Color',[1, 1 ,1])
% 
% % Probability of detection plot
% figure
% imagesc([emitter_bounds(1) emitter_bounds(2)],[emitter_bounds(3) ...
%     emitter_bounds(4)], prob_detection); hold all
% h = colorbar;           % displays the colorbar legend
% colormap('jet')         % sets the colors used
% caxis([0 1])           % sets the limits of the colormap
% set(gca,'YDir','normal')% keeps y-axis correct orientation
% xlabel('x (m)')
% ylabel('y (m)')
% title('Probability of Detection')
% title(h,sprintf('Probability'))
% 
% for ii = 1:numrefs
%     plot(refPos(1,ii), refPos(2,ii), 'k.','HandleVisibility','off', ...
%         'MarkerFaceColor', [0 1 0], 'MarkerSize',28);
%     h = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
%         'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
%         'FontSize', 6);
%     set(h, 'Color',[1, 1 ,1])
% end
% h = text(emitter_bounds(1),0.9*emitter_bounds(4), ...
%     sprintf('Tx EIRP: %2.0f dBm', tx_pwr_dbm));
% set(h, 'Color',[1, 1 ,1])
% 
% % Probability of correlation plot
% figure
% imagesc([emitter_bounds(1) emitter_bounds(2)],[emitter_bounds(3) ...
%     emitter_bounds(4)], prob_correlation); hold all
% h = colorbar;           % displays the colorbar legend
% colormap('jet')         % sets the colors used
% caxis([0 1])           % sets the limits of the colormap
% set(gca,'YDir','normal')% keeps y-axis correct orientation
% xlabel('x (m)')
% ylabel('y (m)')
% title('Probability of Correlation')
% title(h,sprintf('Probability'))
% 
% for ii = 1:numrefs
%     plot(refPos(1,ii), refPos(2,ii), 'k.','HandleVisibility','off', ...
%         'MarkerFaceColor', [0 1 0], 'MarkerSize',28);
%     h = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
%         'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
%         'FontSize', 6);
%     set(h, 'Color',[1, 1 ,1])
% end
% h = text(emitter_bounds(1),0.9*emitter_bounds(4), ...
%     sprintf('Tx EIRP: %2.0f dBm', tx_pwr_dbm));
% set(h, 'Color',[1, 1 ,1])
% 
% fprintf(1,'\nTotal Runtime: %2.1f min\n\n', ...
%     run_time);

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
% for ii = 1:numrefs
%     plot(refPos(1,ii), refPos(2,ii), 'k.','HandleVisibility','off', ...
%         'MarkerFaceColor', [0 1 0], 'MarkerSize',28);
%     h = text(refPos(1,ii), refPos(2,ii)+3, sprintf('%i', ii), ...
%         'horizontalalignment', 'center', 'verticalalignment', 'top');
%     set(h, 'Color',[1, 1 ,1])
% end


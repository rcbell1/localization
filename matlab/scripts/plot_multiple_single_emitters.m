clear; close all
addpath('../functions')

plot_toa_countours = 0;
apply_calibration = 0; % for hardware sims
calib_path = [];
Ntrials = 1000;            % keep at 1 for this plot, only for heat maps
% plot_hull = 1;          % plot convex hull or ref stations
% emitter_bounds = [-100 100 -100 100];   % bounds of emitter locations
% plot_bounds = [-250 250 -250 250];  % outer boundary of plot figure

%% Channel parameters
delay_spread = 300e-9;   % time difference between first received path and last (s)
multi_option = 0;   % determines the type of plot to generate
                    % 0 = no multipath
                    % 1 = two path with various delays between them
                    % 2 = varrying number of paths between 1 and num_paths
                    % with increasing delay spread
num_paths = 2;      % number of multipaths per reciever for option 1
multi_jump = 5;     % skip amount for option 1, 1:multi_jump:num_paths
num_delay_spreads = 10; % number of delay spreads to test in option 2
max_num_paths = inf; % max number paths for option 2

% Generate target emitter positions [x; y] (meters)
% numEmitters = 100;       % total number of random emitter locations to create
% targetPos = [randi([emitter_bounds(1) emitter_bounds(2)],1,numEmitters); ...    
%     randi([emitter_bounds(3) emitter_bounds(4)],1,numEmitters)]; 
% [Tx, Ty] = meshgrid(emitter_bounds(1):emitter_spacing:emitter_bounds(2), ...
%     emitter_bounds(3):emitter_spacing:emitter_bounds(4));

%% Needed by dpd
% -43,3,-4,37
adder = 100;
grid_xmin = -43 - adder;
grid_xmax = 3 + adder;
grid_ymin = -4 - adder;
grid_ymax = 37 + adder;
grid_numx = 20;
grid_numy = 20;

grid_def = [grid_xmin grid_xmax;
            grid_ymin grid_ymax;
            grid_numx grid_numy];

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
targetPos = [targetPos1 targetPos2 targetPos3 targetPos4];
[numdims, numtargets] = size(targetPos);
% targetPos = [targetPos1];

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
tx_pwr_dbm = -37;         % emitter transmit power in dBm (USRP max is 10 dBm)
fs_tx = 200e6/10;
Nsym = 1000;              % number of symbols in signals
span = 10;              % total length of shaping filter in symbols
sps = 2;                % samples per symbol at the transmitter
fsym = fs_tx/sps;             % symbol rate of transmitter (signal bandwidth)
beta = 0.4;             % excess bandwidth of tx pulse shaping filter
fc = 2.395e9;             % center frequency of transmitter

%% Receiver properties
fs = 200e6/1;                % receiver sample rates (Hz)
wlen = 2*ceil(fs/fsym)+1; % moving maximum window length in samples, odd number
nstds = 9;                % number of standard deviations to declare peak
percent_of_peak = 0.8;    % get the number of samples needed on either side 
                          % of correlation peaks for the peak value to drop 
                          % by this percent for use in super resolution

% General simulation properties
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
else
    fprintf(1,'\n\nOption not implemented\n\n')
end
num_multi_idxs = length(multi_idxs);
initial_coords = [-41;1];

jj = 1; % for future use
for ii = 1:size(targetPos,2)
    [coords{ii,jj}, bias_coords{ii,jj}, covar_coords{ii,jj}, ...
    mse_coords(ii,jj), tdoas_true(ii,:), tdoas_coarse{ii,jj}, ...
    tdoas_refined{ii,jj}, prob_correlation(ii,jj), prob_detection(ii,jj), ... 
    dpd_grid, a, unique(ii,jj), corr_mag_sq] = get_single_emitter2(targetPos(:,ii), refPos, ...
    Ntrials, tx_pwr_dbm, fc, fs, fsym, Nsym, span, sps, beta, wlen, ...
    nstds, percent_of_peak, apply_calibration, calib_path, grid_def, ...
    delay_spread, num_paths, max_num_paths, multi_idxs, multi_option, ...
    initial_coords, 0);
end
toc

%% Plots
% Plot localization results
figure
numtargets = size(targetPos,2);
for kk = 1:numtargets
    subplot(1,numtargets,kk)
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
    plot(targetPos(1,kk),targetPos(2,kk), 'kx', 'MarkerSize',18);

%     for ii = 1:grid_numx*grid_numy
%         plot(dpd_grid(1,ii), dpd_grid(2,ii), 'k.', 'MarkerFaceColor', 'k', ...
%             'MarkerSize',6, 'HandleVisibility','off'); 
%         hold all
%     end
    
    for nn = 1:length(coords{kk,1})
        plot(coords{kk,1}(1,nn),coords{kk,1}(2,nn), 'b.', 'MarkerSize',12, 'HandleVisibility','off');
    end
    plot(coords{kk,1}(1),coords{kk,1}(2), 'b.', 'MarkerSize',12)
    axis equal
    axis(bounds(kk,:))
    xlabel('x (m)')
    ylabel('y (m)')
    xlen = bounds(kk,2) - bounds(kk,1);
    ylen = bounds(kk,4) - bounds(kk,3);
    text(bounds(kk,2)-0.35*xlen,bounds(kk,4)-0.15*ylen, ...
     sprintf('b_x: %3.2f (m)\nb_y: %3.2f (m)\n\\sigma_x: %3.2e (m)\n\\sigma_y: %3.2e (m)', ...
     bias_coords{kk}(1), bias_coords{kk}(2), sqrt(covar_coords{kk}(1,1)), ...
     sqrt(covar_coords{kk}(2,2))), 'fontsize', 8)
%     title('Localization Results')
    
    % Plot time contours
    if plot_toa_countours == 1
        tbounds = 2*bounds(kk,:);
        steps = .006*max(tbounds);
        [tx,ty] = meshgrid(tbounds(1):steps:tbounds(2),tbounds(3):steps:tbounds(4));
        tx2 = tx - targetPos(1,kk);
        ty2 = ty - targetPos(2,kk);
        tz = sqrt(tx2.^2+ty2.^2)/c*1e9;
        contour(tx,ty,tz, 4:15, '--', 'showtext', 'on', 'HandleVisibility', 'off')
    end
end
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'Position', [0.5 0.9 0.1 0.1])

figure
xmin = inf;
xmax = -inf;
for nn = 1:numrefs-1
    refined_tdoas = tdoas_refined{kk,1}(:,nn)*1e9;
    mean_tdoa = mean(refined_tdoas);
    std_tdoa = std(refined_tdoas);
    htrue = plot(tdoas_true(kk,nn), nn, 'kx', 'markersize', 8); hold all
    hcoarse = plot(tdoas_coarse{kk,1}(:,nn)*1e9, nn, 'ro');
    hrefine = plot(tdoas_refined{kk,1}(:,nn)*1e9, nn, 'b.');
    xmin = min([xmin; tdoas_true(kk,nn); tdoas_refined{kk,1}(:,nn)*1e9]);
    xmax = max([xmax; tdoas_true(kk,nn); tdoas_refined{kk,1}(:,nn)*1e9]);
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

%% Create table of TDOAs
% for kk = 1:numtargets
%     figure
%     Coarse = tdoas_coarse{kk,1}*1e9;
%     Fine = tdoas_refined{kk,1}*1e9;
%     T1 = table(Coarse, Fine);
%     uitable('ColumnWidth',{108 108 108 108}, 'Data',T1{:,:},'ColumnName',{'Coarse 12','Coarse 13','Fine 12','Fine 13'},...
%     'RowName',T1.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 0.85]);
% 
%     True = tdoas_true(kk,:);
%     T2 = table(True);
%     uitable('ColumnWidth',{108 108}, 'Data',T2{:,:},'ColumnName',{'True 12','True 13'},...
%     'RowName',T2.Properties.RowNames,'Units', 'Normalized', 'Position',[0.2, 0.85, 0.4, 0.15]);
%     % annotation('textbox', [0.2, 0.9, 0.1, 0.1], 'String', sprintf('%3.2f\n',tdoas_true));
% end


% figure
% numrefs = size(refPos,2);
% for ii = 1:numrefs
%     plot(refPos(1,ii), refPos(2,ii), 'r*','HandleVisibility','off', 'MarkerSize',8); 
%     hold all
%     text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
% end
% plot(refPos(1,1), refPos(2,1), 'r*')
% 
% for ii = 1:size(coords,2)
%     plot(targetPos(1,ii),targetPos(2,ii), 'bo', 'MarkerSize',8, 'HandleVisibility','off');
%     plot(coords(1,ii),coords(2,ii), 'b.', 'MarkerSize',12, 'HandleVisibility','off');
% end
% plot(targetPos(1,1),targetPos(2,1), 'bo')
% plot(coords(1,1),coords(2,1), 'b.', 'MarkerSize',12)
% axis equal
% axis(plot_bounds)
% xlabel('x (m)')
% ylabel('y (m)')
% title('Localization Results')
% legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'location', 'best')
% grid on
% 
% if plot_hull == 1
%     k = convhull(refPos(1,:),refPos(2,:));
% end
% plot(refPos(1,k),refPos(2,k),'r', 'HandleVisibility','off')

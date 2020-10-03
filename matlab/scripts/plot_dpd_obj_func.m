clear; close all
addpath('../functions')

%% Sim parameters  
Ntrials = 10;
show_plots = 0;
apply_calibration = 0; % for hardware sims
calib_path = [];

%% Channel parameters
delay_spread = 300e-9;   % time difference between first received path and last (s)
multi_option = 0;   % determines the type of plot to generate
                    % 0 = no multipath
                    % 1 = two path with various delays between them
                    % 2 = varrying number of paths between 1 and num_paths
                    % with increasing delay spread
num_paths = 2;      % number of multipaths per reciever for option 1
multi_jump = 4;     % skip amount for option 1, 1:multi_jump:num_paths
num_delay_spreads = 10; % number of delay spreads to test in option 2
max_num_paths = inf; % max number paths for option 2
max_nlos_amp = 2;
min_num_taps = 100;
multi_dist_based = 0;   % is the NLOS amp based on distance traveled?

%% Emitter pulse properties
tx_pwr_dbm = 67;         % emitter transmit power in dBm (USRP max is 10 dBm)
fs_tx = 200e6/2;
Nsym = 10;              % number of symbols in signals
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
                          
%% Receiver coords using polar coords
a = 39.4908;
r1 = [22.8, -30*pi/180];
r2 = [22.8, -150*pi/180];
r3 = [22.8, 90*pi/180];
center = r1;
refPos = [ r1(1)*cos(r1(2)) r2(1)*cos(r2(2)) r3(1)*cos(r3(2)); ...   % equilateral triangle
           r1(1)*sin(r1(2)) r2(1)*sin(r2(2)) r3(1)*sin(r3(2))];
center = [center(1)*cos(center(2));center(1)*sin(center(2))];
center = [0;0];
refPos = -center + refPos; % origin centered equilateral triangle

%% Emiiter coords
targetPos1 = [sum(refPos(1,:))/3; sum(refPos(2,:))/3];     % center
targetPos2 = [-a/4; sqrt(3)*a/4];     % east
targetPos3 = [-3*a/4; sqrt(3)*a/4];   % west
targetPos4 = [-a/2; 0];               % south
targetPos = [targetPos1];
targetPos = [-1;0];
p0 = targetPos + 0.1*randn(2,1);
[Ndims, Ntx] = size(targetPos);

%% The figures will be bounded by this region
bcenter = [sum(refPos(1,:))/3; sum(refPos(1,:))/3; sum(refPos(2,:))/3; sum(refPos(2,:))/3].';
% bounds = bcenter + [-3 3 -3 3;
%                   -3 3 -3 3;
%                   -3 3 -3 3
%                   -3 3 -3 3];

temp = ones(2,2*Ntx);
temp(1,:) = -temp(1,:);
temp2 = reshape(temp,4,[]).';
r1 = 25;
a = sqrt(2*r1(1)^2+4*r1(1)*cos(120*pi/180));
bounds = bcenter + 1.0*a*reshape(temp,4,[]).';

%% Needed by dpd, defines 2D grid to search over
xmaxref = 1.2*max(refPos(1,:));
xminref = 1.2*min(refPos(1,:));
ymaxref = 1.2*max(refPos(2,:));
yminref = 1.2*min(refPos(2,:));
grid_bounds = [xminref xmaxref yminref ymaxref];
adder = 0;
grid_numx = 30;
grid_numy = 30;
grid_xmin = grid_bounds(1) + bcenter(1) - adder;
grid_xmax = grid_bounds(2) + bcenter(1) + adder;
grid_ymin = grid_bounds(3) + bcenter(3) - adder;
grid_ymax = grid_bounds(4) + bcenter(3) + adder;

grid_def = [grid_xmin grid_xmax;
            grid_ymin grid_ymax;
            grid_numx grid_numy];


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

tic
%% Start processing
% Generate the signal emitted by the target
[x, noise_bw] = generate_signal2(Nsym, fsym, sps, span, beta, show_plots);

% Resample from tx sample rate to rx sample rate
[P,Q] = rat(fs/(fsym*sps));
y1 = resample(x, P, Q);

% Add proper delays that correspond to target and emitter locations
[y2, toas_true, tdoas_true, ranges] = add_delay2(y1, targetPos, refPos, ...
    fs, show_plots);

for nn = 1:Ntrials
    % Add multipath
    [y3, a(:,nn)] = add_multipath(y2, fc, fs, ranges, delay_spread, num_paths, ...
        max_num_paths, multi_idxs{1}, multi_option, max_nlos_amp, ...
        min_num_taps, multi_dist_based, multi_delays, show_plots);
    
    % Add noise at the proper SNR levels for free space path losses
    y4 = add_noise(y3, tx_pwr_dbm, noise_bw, fc, ranges, show_plots);
%     [toas, tdoas] = get_true_toas(refPos, targetPos4);
%     toas
%     [coords(:,nn), dpd_grid, objective_vals{nn}, unique] = dpd(y4, fs, refPos, grid_def);
    [coords(:,nn), dpd_grid, objective_vals{nn}, unique] = ...
        dpd_direct(y4, fs, fc, refPos, grid_def, y1, tx_pwr_dbm, p0);
    
end
toc

[xidx, yidx] = find(objective_vals{1} == min(min(objective_vals{1})));

%% Plot localization results
figure
for kk = 1:num_kk
    % plot the rx locations
    numrefs = size(refPos,2);
    for ii = 1:numrefs
        legh(1) = plot(refPos(1,ii), refPos(2,ii), 'ks', 'MarkerFaceColor', 'k', ...
            'MarkerSize',10, 'HandleVisibility','off'); 
        hold all
    end
    
    % add rx labels
    for ii = 1:numrefs
        h = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle',...
            'FontSize', 8);
        set(h, 'Color',[1, 1 ,1])
    end
    
    % plot the target location
    legh(2) = plot(targetPos(1),targetPos(2), 'kx', 'MarkerSize',12);
    
    % plot estimated emitter location
    for ii = 1:Ntrials
        legh(3) = plot(coords(1,ii),coords(2,ii), 'b.', 'MarkerSize',...
            6, 'HandleVisibility','off');
    end
    axis(bounds(1,:))
    axis equal
    xlabel('x (m)')
    ylabel('y (m)')
    
    % plot the initial estimate
    legh(4) = plot(p0(1), p0(2), 'o');
    
    %     % plot dpd grid
%     dpdh = plot(dpd_grid{1}, dpd_grid{2}, 'k.', 'MarkerFaceColor', 'k', ...
%             'MarkerSize',6, 'HandleVisibility','off'); 
%     legh(4) = dpdh(1);
end
% legend(legh, 'Receiver Locations', 'Target Emitter Location', 'DPD Search Grid', 'Estimated Target Location', 'Location', 'South')
legend(legh, 'Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'Initial Estimate', 'Location', 'South')
% legend(legh, 'Receiver Locations', 'Target Emitter Location', 'DPD Search Grid', 'Estimated Target Location', 'Position', [0.5 0.9 0.1 0.1])

%% Plot objective function heatmap
figure
for kk = 1:num_kk
    % plot objective function
    hpc = pcolor(dpd_grid{1}, dpd_grid{2}, objective_vals{1}); hold all
    colormap('jet')
    % shading interp;
    h = colorbar;
    caxis('auto')
    set(gca,'YDir','normal','color', 'w')% keeps y-axis correct orientation
    set(hpc, 'EdgeColor', 'none')
    % hpc.FaceColor = 'interp';
    xlabel('x (m)')
    ylabel('y (m)')
    axis tight
    title('Objective Function Heatmap')
    title(h, sprintf('Obj Val'))
    
    % plot the rx locations
    numrefs = size(refPos,2);
    for ii = 1:numrefs
        plot(refPos(1,ii), refPos(2,ii), 'ks', 'MarkerFaceColor', 'k', ...
            'MarkerSize',10, 'HandleVisibility','off'); 
    end
    
    % add rx labels
    for ii = 1:numrefs
        h = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle',...
            'FontSize', 8);
        set(h, 'Color',[1, 1 ,1])
    end
    
    % plot the target location
    plot(targetPos(1),targetPos(2), 'wo', 'MarkerSize',6);

    % plot dpd grid
    dpdh = plot(dpd_grid{1}, dpd_grid{2}, 'w.', 'MarkerFaceColor', 'k', ...
            'MarkerSize',6, 'HandleVisibility','off'); 
    legh(3) = dpdh(1);
    
    % plot minimum obj val location
    dpdh = plot(dpd_grid{1}(xidx,yidx), dpd_grid{2}(xidx,yidx), 'wx', ...
            'MarkerSize',6, 'HandleVisibility','off'); 
end
% legend(legh, 'Receiver Locations', 'Target Emitter Location', 'DPD Search Grid', 'Estimated Target Location', 'Location', 'South')

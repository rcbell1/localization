% the number of paths includes the line-of-sight (LOS) path in this script
% for example if num_paths = 2, then there is the direct path and one
% non-LOS (NLOS) path
clear; close all
addpath('../functions')

%% General simulation properties
% apply_calibration = 0; % for hardware sims
% calib_path = [];
Ntrials = 10;            
plot_toa_contours = 0;                      % 0 off, 1 on
sim_params.Ntrials = Ntrials;
sim_params.apply_calibration = 0;           % for hardware sims
sim_params.calib_path = [];                 % for hardware sims
sim_params.show_plots = 0;
loc_types = {'sparse-dpd','dpd','lsq', 'si', 'taylor'};
% loc_types = {'dpd','lsq', 'si', 'taylor'};

%% Channel parameters
delay_spread = 300e-8;   % time difference between first received path and last (s)
multi_option = 2;   % determines the type of plot to generate
                    % 0 = no multipath
                    % 1 = two path with various delays between them
                    % 2 = varrying number of paths between 1 and max_num_paths
                    % with increasing delay spread
                    % 3 = two path with increasing delay spread
num_paths = 2;      % number of multipaths per reciever for option 1
multi_jump = 4;     % skip amount for option 1, 1:multi_jump:num_paths
num_delay_spreads = 10; % number of delay spreads to test in option 2
max_num_paths = inf; % max number paths for option 2
max_nlos_amp = 1;
min_num_taps = 100;
multi_dist_based = 1;   % is the NLOS amp based on distance traveled?

channel_params.multi_option = multi_option;
channel_params.num_paths = num_paths;
channel_params.multi_jump = multi_jump;
channel_params.num_delay_spreads = num_delay_spreads;
channel_params.max_num_paths = max_num_paths;
channel_params.max_amp = max_nlos_amp;
channel_params.min_num_taps = min_num_taps;
channel_params.distance_based = multi_dist_based;

%% Receiver coords by spacing rxs uniformly on circle
Nrx = 4;
radius = 80;
origin = 0;     % receiver number to be assigned coordinate [0;0], 0 means none
[refPos, center] = get_rx_coords(Nrx, origin, radius);
receiver_params.ref_locs = refPos;
multi_coords = repmat({nan},Nrx,1);
multi_coords{1} = [-5; 30];
% multi_coords{1} = [0; 0];
% multi_coords{2} = [-80; 80];
% multi_coords{3} = [-40 0; 40 5];
% multi_coords{2} = [-20; 0];
channel_params.multi_coords = multi_coords;

%% Emitter coords rectangular coords
targetPos1 = [25.15;-50.15];
targetPos1 = [3;27];
targetPos1 = [0;0];
% targetPos1 = [randi([-50 50]);randi([-50 50])];
targetPos = [targetPos1] + center;

[transmitter_params.Ndims, transmitter_params.Ntargets] = size(targetPos);
transmitter_params.target_locs = targetPos;

%% The figures will be bounded by this region
% Nemitter = size(targetPos,2);
fig_bounds = ones(2,2);
fig_bounds(1,:) = -fig_bounds(1,:);
fig_bounds = center' + 1.5*radius*fig_bounds;
              
%% Emitter pulse properties
tx_pwr_dbm = -10:2:20;         % emitter transmit power in dBm (USRP max is 10 dBm)
% fs_tx = 200e6/2.5; %5
fs_tx = 200e6/20;
Nsym = 10;              % number of symbols in signals
span = 10;              % total length of shaping filter in symbols
sps = 2;                % samples per symbol at the transmitter
fsym = fs_tx/sps;             % symbol rate of transmitter (signal bandwidth)
Tsym = sps/fs_tx;
beta = 0.4;             % excess bandwidth of tx pulse shaping filter
fc = 2.395e9;             % center frequency of transmitter

transmitter_params.sample_rate = fs_tx;
transmitter_params.Nsyms = Nsym;
transmitter_params.span = span;
transmitter_params.sps = sps;
transmitter_params.symbol_rate = fsym;
transmitter_params.symbol_period = Tsym;
transmitter_params.excess_bw = beta;
transmitter_params.carrier_freq = fc;

%% Receiver properties
fs = 200e6/10;                % receiver sample rates (Hz)
percent_of_peak = 0.8;    % get the number of samples needed on either side 
                          % of correlation peaks for the peak value to drop 
                          % by this percent for use in super resolution
                          
receiver_params.sample_rate = fs;     % receiver sample rates (Hz)
receiver_params.percent_of_peak = percent_of_peak; % get the number of samples 
                            % needed on either side of correlation peaks 
                            % for the peak value to drop by this percent 
                            % for use in super resolution
%% Needed by dpd
% brad = 1.2*radius;
brad = 1*radius;
grid_bounds = [-brad brad -brad brad];
adder = 0;
grid_center = targetPos;
grid_center = [0;0];
grid_xmin = grid_bounds(1) - adder + grid_center(1);
grid_xmax = grid_bounds(2) + adder + grid_center(1);
grid_ymin = grid_bounds(3) - adder + grid_center(2);
grid_ymax = grid_bounds(4) + adder + grid_center(2);
grid_numx = 31;
grid_numy = 31;

grid_def = [grid_xmin grid_xmax;
            grid_ymin grid_ymax;
            grid_numx grid_numy];

loc_alg_params.grid_def = grid_def;

start_time = tic;
tic
num_ds_samps = ceil(delay_spread*fs); % delay spread in rx samples count
if multi_option == 0
    delay_spread = zeros(num_delay_spreads,1);
    multi_idxs = cell(num_ds_samps,1);
    num_jj = 1;
elseif multi_option == 1
    delay_spread = repmat(delay_spread,num_delay_spreads,1);
    multi_idxs = num2cell(0:multi_jump:num_ds_samps);
    num_jj = length(multi_idxs);
elseif multi_option == 2
    delay_spread = linspace(0,delay_spread,num_delay_spreads);
    multi_idxs = cell(num_delay_spreads,1); % this will be filled in randomly later
    multi_idxs(:) = {nan};
    num_jj = num_delay_spreads;
elseif multi_option == 3
    delay_spread = linspace(0,delay_spread,num_delay_spreads);
    multi_idxs = cell(num_delay_spreads,1); % this will be filled in randomly later
    multi_idxs(:) = {nan};
    num_jj = num_delay_spreads;
elseif multi_option == 4
    delay_spread = zeros(num_delay_spreads,1);
    multi_idxs = cell(num_ds_samps,1);
    num_jj = 1;
else
    fprintf(1,'\n\nOption not implemented\n\n')
end
initial_coords = targetPos + [-10;10] + center;
loc_alg_params.initial_coords = initial_coords;
Ntypes = length(loc_types);
Npowers = length(tx_pwr_dbm);
for ii = 1:Ntypes
    loc_alg_params.type = loc_types{ii};
    for jj = 1:num_jj 
        channel_params.multi_idxs = multi_idxs{jj};
        num_multi_idxs = length(multi_idxs);
        for kk = 1:Npowers
        
            channel_params.delay_spread = delay_spread(jj);
            transmitter_params.tx_pwr_dbm = tx_pwr_dbm(kk);

            [results.target{ii}.multi{jj}.power{kk}] = ...
            get_single_emitter3(sim_params, transmitter_params, ...
            receiver_params, channel_params, loc_alg_params);

            % Unpack results to local variables
            coords{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.coords;
            bias_coords{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.bias_coords;
            covar_coords{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.covar_coords;
            mse_coords(ii,jj,kk) = results.target{ii}.multi{jj}.power{kk}.mse_coords;
            avg_snr_db{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.avg_snr_db;
            h{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.channel_taps;
            tx_raw{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.tx_raw;
            tx_up{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.tx_up;
            tx_delayed{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.tx_delayed;
            tx_multi{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.tx_multi;
            tx_noise{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.tx_noise;
            tx_agc{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.tx_agc;
            multi_delays = results.target{ii}.multi{jj}.power{kk}.multi_delays;
            ranges = results.target{ii}.multi{jj}.power{kk}.ranges;
            tdoas_true(:,ii,kk) = results.target{ii}.multi{jj}.power{kk}.tdoas_true;
            toas_true(:,ii,kk) = results.target{ii}.multi{jj}.power{kk}.toas_true;

            switch loc_alg_params.type
                case {'dpd','direct-dpd','ms-dpd1','ms-dpd2','ms-music-dpd','ms-music-dpd2','sparse-dpd'}
                    if contains(loc_alg_params.type, 'ms-dpd')
                        rcoords{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.reflection_coords;
                    end
                    dpd_grid{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.grid;
                    if ~contains(loc_alg_params.type, 'sparse-dpd')
                        dpd_obj_vals{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.obj_vals;
                    end
                    if contains(loc_alg_params.type, 'sparse-dpd')
                        sparse_heatmap{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.sparse_heatmap;
                        sparse_max_coord{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.max_coord;
                        sparse_all_coords{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.all_coords;
                    end

                case {'lsq', 'si', 'taylor'}
                    tdoas_coarse{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.tdoas_coarse;
                    tdoas_refined{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.tdoas_refined;
                    prob_correlation(ii,jj,kk) = results.target{ii}.multi{jj}.power{kk}.prob_correlation;
                    prob_detection(ii,jj,kk) = results.target{ii}.multi{jj}.power{kk}.prob_detection;
                    corr_mag{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.corr_mag;
                    corr_mag_bad{ii,jj,kk} = results.target{ii}.multi{jj}.power{kk}.corr_mag_bad;

                    if strcmp(loc_alg_params.type, 'lsq')
                        unique(ii,jj,kk) = results.target{ii}.multi{jj}.power{kk}.unique;                    
                    end
            end % switch
        end
    end
end
toc
stop_time = toc(start_time);

%% New Plots
markers = {'x';'o';'d';'+';'s';'.'};
leg_pos = [0.2 0.67 0.1 0.1]; % legend position
fig_dims = [300, 35, 1100, 725]; % figure window defs
anno_pos = [0.1 0.75 0.2 0.2]; % annotation box position
my_anno = [];
for ii = 1:max(Nrx,6)
    if ii > Nrx
        my_anno{ii} = sprintf('                ');
    else        
        my_anno{ii} = sprintf('SNR%i:%3.1f', ii, avg_snr_db{1,1}(ii,1));
    end
end
my_anno{1} = [my_anno{1} ', DelSprd:' ...
    sprintf('%3.0f/%i',delay_spread(end)*1e9,num_delay_spreads) ...
    ', MDistBased:' sprintf('%i',multi_dist_based)];
my_anno{2} = [my_anno{2} ', Tx BW:' sprintf('%2.1f',fsym/1e6) ', Ntrials:' sprintf('%i',Ntrials)];
my_anno{3} = [my_anno{3} ', Rx Samp:' sprintf('%2.1f',fs/1e6) ', MOption:' sprintf('%i',multi_option)];
my_anno{4} = [my_anno{4} ', MPaths:' sprintf('%i',num_paths) ', Nsym:' sprintf('%i',Nsym)] ;
my_anno{5} = [my_anno{5} ', MaxMAmp:' sprintf('%i',max_nlos_amp) ', fc:' sprintf('%3.3f',fc*1e-9)];
my_anno{6} = [my_anno{6} ', simtime:' sprintf('%.1f',stop_time/60)];

figure
annotation('textbox',anno_pos,'String',my_anno, 'FitBoxToText','on');
for jj = 1:num_jj
    if num_jj > 1
        subplot(ceil((num_jj+2)/4),4,jj+2)
    else
        subplot(ceil((num_jj+2)/4),1,jj+2)
    end
    mse_plot = squeeze(mse_coords(:,jj,:));
    hf = plot(tx_pwr_dbm, 10*log10(mse_plot));
    hf = plot(tx_pwr_dbm, mse_plot);
    set(hf,{'Marker'},markers(1:Ntypes))
%     axis([-inf inf 0 50])
    axis([-inf inf -inf inf])
    if num_jj > 1
        title( sprintf('%4.1f',delay_spread(jj)*1e9) )
    end
    legend(hf, loc_types{:}, 'Position', leg_pos)
end
set(gcf, 'Position',  fig_dims)

% Plot localization results
figure
if num_jj == 1
    blanks = 1;
    if contains([loc_types{:}], 'sparse-dpd')
        num_subplots = 3; 
        max_cols = 3;
    else
        num_subplots = 2;
        max_cols = 2;
    end
    num_rows = 1;
    leg_pos = [0.25 0.4 0.1 0.1]; % legend position
    fig_dims = [300, 100, 1000, 500]; % figure window defs
    anno_pos = [0.1 0.65 0.2 0.2]; % annotation box position
    
%     subplot(num_rows,max_cols,kk)
    numrefs = size(refPos,2);
    for ii = 1:numrefs
        tmp = plot(refPos(1,ii), refPos(2,ii), 'ks', 'MarkerFaceColor', 'k', ...
            'MarkerSize',10, 'HandleVisibility','off'); 
        hold all
    end
    lh = tmp(1);
    
    % plotting text around markers
    for ii = 1:numrefs
        ht = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle',...
            'FontSize', 8);
        set(ht, 'Color',[1, 1 ,1])   
    end
    
    % Plot the true emitter position
    tmp = plot(targetPos(1),targetPos(2), 'bp', 'MarkerSize',10,'MarkerFaceColor', [70,190,238]/255);
    lh = [lh tmp(1)];

    % plot dpd grid
    if contains([loc_types{:}], 'dpd')
        tmp = plot(dpd_grid{1,1}{1}, dpd_grid{1,1}{2}, 'k.', 'MarkerFaceColor', 'k', ...
                'MarkerSize',6, 'HandleVisibility','off'); 
        lh = [lh tmp(1)];
    end
    
    if contains([loc_types{:}], 'taylor')
        tmp = plot(initial_coords(1), initial_coords(2), 'ko');
        lh = [lh tmp(1)];
    end
    
    axis equal
    axis(fig_bounds(:))
    if strcmp(loc_alg_params.type,'dpd')
        axis tight
    end

    xlabel('x (m)')
    ylabel('y (m)')
    set(gcf, 'Position',  fig_dims)
end


%% Old Plots
c = 299792458;          % speed of light m/s
my_anno = [];
for ii = 1:max(Nrx,6)
    if ii > Nrx
        my_anno{ii} = sprintf('                ');
    else        
        my_anno{ii} = sprintf('SNR%i:%3.1f', ii, avg_snr_db{1,1}(ii,1));
    end
end
my_anno{1} = [my_anno{1} ', Alg:' loc_alg_params.type ', DelSprd:' ...
    sprintf('%3.0f/%i',delay_spread(end)*1e9,num_delay_spreads) ...
    ', MDistBased:' sprintf('%i',multi_dist_based)];
my_anno{2} = [my_anno{2} ', Tx BW:' sprintf('%2.1f',fsym/1e6) ', Ntrials:' sprintf('%i',Ntrials)];
my_anno{3} = [my_anno{3} ', Rx Samp:' sprintf('%2.1f',fs/1e6) ', MOption:' sprintf('%i',multi_option)];
my_anno{4} = [my_anno{4} ', MPaths:' sprintf('%i',num_paths) ', Nsym:' sprintf('%i',Nsym)] ;
my_anno{5} = [my_anno{5} ', MaxMAmp:' sprintf('%i',max_nlos_amp) ', fc:' sprintf('%3.3f',fc*1e-9)];
my_anno{6} = [my_anno{6} ', TxPdBm:' sprintf('%2.1f',tx_pwr_dbm) ', simtime:' sprintf('%.1f',stop_time/60)];
if multi_option == 1
    
    temp = horzcat(bias_coords{:});
    bx = temp(1,:);
    by = temp(2,:);
    temp = horzcat(covar_coords{:});
    for jj = 1:num_multi_idxs
        temp = diag(covar_coords{1,jj});
        stdx(jj) = sqrt(temp(1));
        stdy(jj) = sqrt(temp(2));
    end

    figure
    annotation('textbox',[0.35 0.75 0.2 0.2],'String',my_anno, 'FitBoxToText','on');
    
    % Plot multipath index vs statistics
    subplot(412)
    plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), sqrt(mse_coords(1,:)), '.-'); hold all
    ylims = ylim;
    ylims(2) = min(60,ylims(2));
    plot([Tsym:Tsym:delay_spread; Tsym:Tsym:delay_spread]/1e-9, [ylims(1) ylims(2)], 'k--');
    title(sprintf('Two Paths with Increasing Delay Between Them, Tx Power %3.0f dBm', tx_pwr_dbm))
    xlabel('Delay Between 1st and 2nd Path (ns)')
    ylabel('RMSE (m)')
    grid on
    axis([-inf inf ylims(1) ylims(2)])

    subplot(413)
    plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), bx, '.-'); hold all
    plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), by, '.-')
    ylims = ylim;
    ylims(1) = max(-60,ylims(1));
    ylims(2) = min(60,ylims(2));
    plot([Tsym:Tsym:delay_spread; Tsym:Tsym:delay_spread]/1e-9, [ylims(1) ylims(2)], 'k--');
    xlabel('Delay Between 1st and 2nd Path (ns)')
    ylabel('Bias (m)')
    legend('b_x','b_y')
    grid on
    axis([-inf inf ylims(1) ylims(2)])

    subplot(414)
    plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), stdx, '.-'); hold all
    plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), stdy, '.-')
    ylims = ylim;
    ylims(2) = min(60,ylims(2));
    plot([Tsym:Tsym:delay_spread; Tsym:Tsym:delay_spread]/1e-9, [ylims(1) ylims(2)], 'k--');
    xlabel('Delay Between 1st and 2nd Path (ns)')
    ylabel('Standard Deviation (m)')
    legend('\sigma_x','\sigma_y')
    grid on
    axis([-inf inf ylims(1) ylims(2)])
    set(gcf, 'Position',  [300, 100, 1100, 800])

    % plot all correlation peaks *************************************
    if ~contains(loc_alg_params.type, 'dpd')
        plot_opts = {'b.-', 'r.-', 'g.-', 'm.-'};
        num_opts = length(plot_opts);
        kk = 1;
        figure
        for ii = 1:num_multi_idxs
            subplot(num_multi_idxs,1,ii)
            nsamps = length(corr_mag{ii});
            lags = (1:nsamps) - (nsamps+1)/2;
            plot(lags/(fs*1e-9), corr_mag{ii}, plot_opts{kk}); hold on
            axis([-3*delay_spread(1) 3*delay_spread(1) -inf inf]/1e-9)

            if kk == num_opts
                kk = 1;
            else
                kk = kk + 1;
            end
        end
        set(gcf, 'Position',  [300, 100, 1300, 800])
    end
    
    % plot the bad correlation peaks *************************************
    if ~contains(loc_alg_params.type, 'dpd')
        t_label = 'Rx Pair 1,';
        for ii = 1:num_delay_spreads
            figure
            %             rx_pair = 1;
            for jj = 1:Nrx-1 
                subplot(1,Nrx-1,jj)
                base = 0;
                sigs = horzcat(corr_mag_bad{ii}{:,jj});
                for kk = 1:size(sigs,2)
                    sig = sigs(:,kk);
                    nsamps = length(sig);
                    lags = ((1:nsamps) - (nsamps+1)/2).';
                    hc = max(abs(sig));
                    if kk == 1
                        plot(lags/(fs*1e-9), real(sig), '.-'); hold all
                        plot(xlim,[base base],'k-')
                        base = base + 1.0*hc;
                        text(-delay_spread(end)*1e9, base, sprintf('%i',kk))
                    else
                        plot(lags/(fs*1e-9), real(sig) + base + hc, '.-');
                        plot(xlim,[base+hc base+hc],'k-')
                        base = base + 2.0*hc;
                        text(-delay_spread(end)*1e9, base, sprintf('%i',kk))
                    end  
                end
                plot([tdoas_true(jj) tdoas_true(jj)]*1e9,ylim,'k--')
                title(sprintf('%s %i', t_label, jj+1))
                xlabel('Lag (ns)')
                ylabel('Bad XCorr Magnitude')
                axis([-1*delay_spread(end) 1*delay_spread(end) -inf inf]/1e-9)
            end
            annotation('textbox',[0.0 0.80 0.1 0.1], 'String', ...
                sprintf('Delay Spread\n%3.0f', delay_spread(ii)*1e9), 'FitBoxToText','on')
            set(gcf, 'Position',  [200, 100, 1300, 800])
        end
    end

elseif multi_option == 2 || multi_option == 3
    
    temp = horzcat(bias_coords{:});
    bx = temp(1,:);
    by = temp(2,:);
    temp = horzcat(covar_coords{:});
    for jj = 1:num_delay_spreads
        temp = diag(covar_coords{1,jj});
        stdx(jj) = sqrt(temp(1));
        stdy(jj) = sqrt(temp(2));
    end
 
    figure
    annotation('textbox',[0.35 0.75 0.2 0.2],'String',my_anno, 'FitBoxToText','on');
    
    subplot(412)
    plot(delay_spread/1e-9, sqrt(mse_coords(1,:)), '.-'); hold all
    ylims = ylim;
%     ylims(2) = min(60,ylims(2));
    if Tsym < max(delay_spread)
        plot([Tsym:Tsym:max(delay_spread); Tsym:Tsym:max(delay_spread)]/1e-9, [0 1.1*ylims(2)], 'k--');
    end
    title(sprintf('Random Number of Paths and Path Delays, Tx Power %3.0f dBm', tx_pwr_dbm))
    xlabel('Delay Spread (ns)')
    ylabel('RMSE (m)')
    grid on
%     axis([-inf inf ylims(1) ylims(2)])
    axis([-inf inf 0 1.1*ylims(2)])

    subplot(413)
    plot(delay_spread/1e-9, bx, '.-'); hold all
    plot(delay_spread/1e-9, by, '.-')
    ylims = ylim;
    ylims(1) = max(-60,ylims(1));
    ylims(2) = min(60,ylims(2));
    if Tsym < max(delay_spread)
        plot([Tsym:Tsym:max(delay_spread); Tsym:Tsym:max(delay_spread)]/1e-9, [ylims(1) ylims(2)], 'k--');
    end
    xlabel('Delay Spread (ns)')
    ylabel('Bias (m)')
    legend('b_x','b_y')
    grid on
    axis([-inf inf ylims(1) ylims(2)])

    subplot(414)
    plot(delay_spread/1e-9, real(stdx), '.-'); hold all
    plot(delay_spread/1e-9, real(stdy), '.-')
    ylims = ylim;
    ylims(2) = min(60,ylims(2));
    if Tsym < max(delay_spread)
        plot([Tsym:Tsym:max(delay_spread); Tsym:Tsym:max(delay_spread)]/1e-9, [ylims(1) ylims(2)], 'k--');
    end
    xlabel('Delay Spread (ns)')
    ylabel('Standard Deviation (m)')
    legend('\sigma_x','\sigma_y')
    grid on
    axis([-inf inf ylims(1) ylims(2)])
    set(gcf, 'Position',  [300, 100, 1100, 800])

    % plot all correlation peaks *************************************
    if ~contains(loc_alg_params.type, 'dpd')
        plot_opts = {'b.-', 'r.-', 'g.-', 'm.-'};
        num_opts = length(plot_opts);
        kk = 1;
        figure
        for ii = 1:num_delay_spreads
            subplot(num_delay_spreads,1,ii)
            nsamps = length(squeeze(corr_mag{ii}(:,:,1).'));
            lags = (1:nsamps) - (nsamps+1)/2;
            plot(lags/(fs*1e-9), squeeze(corr_mag{ii}(:,:,1).'), plot_opts{kk}); hold on
            axis([-1*delay_spread(end) 1*delay_spread(end) -inf inf]/1e-9)

            if kk == num_opts
                kk = 1;
            else
                kk = kk + 1;
            end
        end
    end
    
    % plot the bad correlation peaks *************************************
    if ~contains(loc_alg_params.type, 'dpd')
        t_label = 'Rx Pair 1,';
        for ii = 1:num_delay_spreads
            figure
            %             rx_pair = 1;
            for jj = 1:Nrx-1 
                subplot(1,Nrx-1,jj)
                base = 0;
                sigs = horzcat(corr_mag_bad{ii}{:,jj});
%                 sig = squeeze(corr_mag_bad{ii}(jj,:,rx_pair).');
                for kk = 1:size(sigs,2)
                    sig = sigs(:,kk);
                    nsamps = length(sig);
                    lags = ((1:nsamps) - (nsamps+1)/2).';
                    hc = max(abs(sig));
                    if kk == 1
                        plot(lags/(fs*1e-9), real(sig), '.-'); hold all
                        plot(xlim,[base base],'k-')
                        base = base + 1.0*hc;
                        text(-delay_spread(end)*1e9, base, sprintf('%i',kk))
                    else
                        plot(lags/(fs*1e-9), real(sig) + base + hc, '.-');
                        plot(xlim,[base+hc base+hc],'k-')
                        base = base + 2.0*hc;
                        text(-delay_spread(end)*1e9, base, sprintf('%i',kk))
                    end  
                end
                plot([tdoas_true(jj) tdoas_true(jj)]*1e9,ylim,'k--')
                title(sprintf('%s %i', t_label, jj+1))
                xlabel('Lag (ns)')
                ylabel('Bad XCorr Magnitude')
                axis([-1*delay_spread(end) 1*delay_spread(end) -inf inf]/1e-9)
            end
            annotation('textbox',[0.0 0.80 0.1 0.1], 'String', ...
                sprintf('Delay Spread\n%3.0f', delay_spread(ii)*1e9), 'FitBoxToText','on')
            set(gcf, 'Position',  [200, 100, 1300, 800])
        end
    end

else
  
    % plot all correlation peaks *************************************
    t_label = 'Rx Pair 1,';
    if ~contains(loc_alg_params.type, 'dpd')
        figure
        for ii = 1:Nrx-1
            subplot(1,Nrx-1,ii)
            nsamps = size(corr_mag{1}(:,:,ii).',1);
            lags = (1:nsamps) - (nsamps+1)/2;
            plot(lags/(fs*1e-9), corr_mag{1}(:,:,ii).', '.-'); hold on
    %         axis([-1000/(fs*1e-9) 1000/(fs*1e-9) -inf inf])
            
            plot([tdoas_true(ii) tdoas_true(ii)]*1e9,ylim,'k--')
            title(sprintf('%s %i', t_label, ii+1))
            xlabel('Lag (ns)')
            ylabel('Bad XCorr Magnitude')
%             axis([-inf inf -inf inf])
            axis([-0.125*nsamps/(fs*1e-9) 0.125*nsamps/(fs*1e-9) -inf inf])
        end
        set(gcf, 'Position',  [200, 100, 1300, 800])

    
        % plot the bad correlation peaks *******************************
        figure
        for jj = 1:Nrx-1 
            subplot(1,Nrx-1,jj)
            base = 0;
            sigs = horzcat(corr_mag_bad{1}{:,jj});
            for kk = 1:size(sigs,2)
                sig = sigs(:,kk);
                nsamps = length(sig);
                lags = ((1:nsamps) - (nsamps+1)/2).';
                hc = max(abs(sig));
                if kk == 1
                    plot(lags/(fs*1e-9), real(sig), '.-'); hold all
                    plot(xlim,[base base],'k-')
                    base = base + 1.0*hc;
                    text(-delay_spread(end)*1e9, base, sprintf('%i',kk))
                else
                    plot(lags/(fs*1e-9), real(sig) + base + hc, '.-');
                    plot(xlim,[base+hc base+hc],'k-')
                    base = base + 2.0*hc;
                    text(-delay_spread(end)*1e9, base, sprintf('%i',kk))
                end  
            end
            plot([tdoas_true(jj) tdoas_true(jj)]*1e9,ylim,'k--')
            title(sprintf('%s %i', t_label, jj+1))
            xlabel('Lag (ns)')
            ylabel('Bad XCorr Magnitude')
%             axis([-1000/(fs*1e-9) 1000/(fs*1e-9) -inf inf])
        end
%         annotation('textbox',[0.0 0.80 0.1 0.1], 'String', ...
%             sprintf('Delay Spread\n%3.0f', delay_spread(1)*1e9), 'FitBoxToText','on')
        set(gcf, 'Position',  [200, 100, 1300, 800])
        axis([-inf inf -inf inf])
    end

end


%% Plot localization results
figure
if num_jj == 1
    blanks = 1;
    if contains(loc_alg_params.type, 'sparse-dpd')
        num_subplots = 3; 
        max_cols = 3;
    else
        num_subplots = 2;
        max_cols = 2;
    end
    num_rows = 1;
    leg_pos = [0.25 0.4 0.1 0.1]; % legend position
    fig_dims = [300, 100, 1000, 500]; % figure window defs
    anno_pos = [0.1 0.65 0.2 0.2]; % annotation box position
else
    blanks = 2; % empty spots for annotation box
    num_subplots = num_jj + blanks; 
    if num_subplots > 4
        max_cols = 4;
    else
        max_cols = num_subplots;
    end
    num_rows = ceil(num_subplots/max_cols);  % 4 plots across max
    leg_pos = [0.15 0.7 0.1 0.1]; % legend position
    fig_dims = [300, 100, 1100, 800]; % figure window defs
    anno_pos = [0.1 0.75 0.2 0.2]; % annotation box position
end
annotation('textbox',anno_pos,'String',my_anno, 'FitBoxToText','on');
for kk = (1+blanks):num_subplots
    subplot(num_rows,max_cols,kk)
    numrefs = size(refPos,2);
    for ii = 1:numrefs
        tmp = plot(refPos(1,ii), refPos(2,ii), 'ks', 'MarkerFaceColor', 'k', ...
            'MarkerSize',10, 'HandleVisibility','off'); 
        hold all
    end
    lh = tmp(1);
    
    % plotting text around markers
    for ii = 1:numrefs
        ht = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle',...
            'FontSize', 8);
        set(ht, 'Color',[1, 1 ,1])   
    end
    
    % '#4DBEEE','#0072BD'
    % Plot the true emitter position
    tmp = plot(targetPos(1),targetPos(2), 'bp', 'MarkerSize',10,'MarkerFaceColor', '#4DBEEE');
    lh = [lh tmp(1)];

    % plot dpd grid
    if contains(loc_alg_params.type, 'dpd')
        tmp = plot(dpd_grid{1,1}{1}, dpd_grid{1,1}{2}, 'k.', 'MarkerFaceColor', 'k', ...
                'MarkerSize',6, 'HandleVisibility','off'); 
        lh = [lh tmp(1)];
    end
    
    if strcmp(loc_alg_params.type, 'sparse-dpd')
        if kk == num_subplots
            for ii = 1:Ntrials
                ecoords = sparse_max_coord{1}(:,ii);
                    tmp = plot(ecoords(1),ecoords(2), 'b.', 'MarkerSize',...
                        14, 'HandleVisibility','off');
            end
        else
            for ii = 1:Ntrials
                ecoords = sparse_all_coords{1}{ii};
                nemitters = size(ecoords,2);
                for jj = 1:nemitters
                    tmp = plot(ecoords(1,jj),ecoords(2,jj), 'b.', 'MarkerSize',...
                        14, 'HandleVisibility','off');
                end
            end
        end
    else
        for ii = 1:Ntrials
            tmp = plot(coords{1,kk-blanks}(1,ii),coords{1,kk-blanks}(2,ii), 'b.', 'MarkerSize',...
                12, 'HandleVisibility','off');
        end
    end
    lh = [lh tmp(1)];
    
    if strcmp(loc_alg_params.type,'taylor')
        tmp = plot(initial_coords(1), initial_coords(2), 'ko');
        lh = [lh tmp(1)];
    end
    
    axis equal
    axis(fig_bounds(:))
    if strcmp(loc_alg_params.type,'dpd')
        axis tight
    end

    title(sprintf('%3.0f ns', delay_spread(kk-blanks)*1e9))
    xlabel('x (m)')
    ylabel('y (m)')
    set(gcf, 'Position',  fig_dims)

    % plot multipath reflection location
    if multi_option == 4
        for ii = 1:length(multi_coords)
            if ~isnan(multi_coords{ii})
                for jj = 1:size(multi_coords{ii},2)
                    coord = multi_coords{ii}(:,jj);
                    plot(coord(1),coord(2), 'k^', 'MarkerSize',10, 'MarkerFaceColor', 'g');
                end
            end
        end

        % add labels
        for ii = 1:length(multi_coords)
            if ~isnan(multi_coords{ii})
                for jj = 1:size(multi_coords{ii},2)
                    coord = multi_coords{ii}(:,jj);
                    h = text(coord(1), coord(2), sprintf('%i', ii), ...
                        'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                        'FontSize', 8);
                    set(h, 'Color',[0, 0 ,0])
                end
            end
        end
    end

    
    % Plot time contours
    if plot_toa_contours == 1
        tbounds = 2*fig_bounds(kk-blanks,:);
        steps = .006*max(tbounds);
        [tx,ty] = meshgrid(tbounds(1):steps:tbounds(2),tbounds(3):steps:tbounds(4));
        tx2 = tx - targetPos(1,kk-blanks);
        ty2 = ty - targetPos(2,kk-blanks);
        tz = sqrt(tx2.^2+ty2.^2)/c*1e9;
        contour(h1,tx,ty,tz, 4:15, '--', 'showtext', 'on', 'HandleVisibility', 'off')
    end
end
switch loc_alg_params.type
    case 'taylor'
        legend(lh, 'Receiver Locations', 'Target Emitter Location', ...
            'Estimated Target Location', 'Initial Guess', 'Position', leg_pos)
    case {'dpd','direct-dpd','ms-dpd1','ms-dpd2','ms-music-dpd','ms-music-dpd2','sparse-dpd'}
        legend(lh, 'Receiver Locations', 'Target Emitter Location', ...
            'DPD Grid Points', 'Estimated Target Location', 'Position', leg_pos)
    otherwise
        legend(lh,'Receiver Locations', 'Target Emitter Location', ...
            'Estimated Target Location', 'Position', leg_pos)
end

%% Plot tdoa estimation and statistics
if ~contains(loc_alg_params.type, 'dpd')
    skips = 2;
    num_subplots = num_jj+skips;
    if num_subplots > 4
        max_rows = 4;
    else
        max_rows = num_subplots;
    end
    num_cols = ceil(num_subplots/max_rows);  % 4 plots across max
    
    figure
    subplot(max_rows,num_cols,2)
    numrefs = size(refPos,2);
    for ii = 1:numrefs
        tmp = plot(refPos(1,ii), refPos(2,ii), 'ks', 'MarkerFaceColor', 'k', ...
            'MarkerSize',10, 'HandleVisibility','off'); 
        hold all
    end
    lh = [];
    lh = tmp(1);
    
    % plotting text around markers
    for ii = 1:numrefs
        ht = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle',...
            'FontSize', 8);
        set(ht, 'Color',[1, 1 ,1]) 
    end
    
    tmp = plot(targetPos(1),targetPos(2), 'kx', 'MarkerSize',18);
    lh = [lh tmp(1)];
    
    axis equal
    axis(fig_bounds(:))
    xlabel('x (m)')
    ylabel('y (m)')

    h2 = gca;
    low_end = max(10,0.3*abs(min(tdoas_true))*1e9);
    high_end = max(10,0.3*abs(max(tdoas_true))*1e9);
    xlims(1) = min(tdoas_true)*1e9 - low_end;
    xlims(2) = max(tdoas_true)*1e9 + high_end;
    for kk = (1+skips):num_subplots
        xmin = inf;
        xmax = -inf;
        xmin = inf;
        xmax = -inf;
        ax(kk-skips) = subplot(max_rows,num_cols,kk);
        mpdelay = multi_idxs{kk-skips}/fs/1e-9;
        for nn = 1:numrefs-1
            refined_tdoas = tdoas_refined{kk-skips}(:,nn)*1e9;
            mean_tdoa = mean(refined_tdoas);
            std_tdoa = std(refined_tdoas);
            htrue = plot(tdoas_true(nn)*1e9, nn, 'kx', 'markersize', 8); hold all
            hcoarse = plot(tdoas_coarse{kk-skips}(:,nn)*1e9, nn, 'ro');
            hrefine = plot(tdoas_refined{kk-skips}(:,nn)*1e9, nn, 'b.');
            text(tdoas_true(nn)*1e9, nn, sprintf('\\mu = %3.2f, \\sigma = %3.2f', mean_tdoa, std_tdoa), ...
                'horizontalalignment', 'center', 'verticalalignment', 'top',...
                'FontSize', 8);
            ylabels{nn} = sprintf('1,%i', nn+1);
        end
        axis([xlims(1) xlims(2) 0.5 numrefs])
        xtick_skip = max(1,floor(length(floor(xlims(1)):floor(xlims(2)))/10));
        xticks(floor(xlims(1)):xtick_skip:floor(xlims(2)))
        y_values = 1:numrefs-1;
        
        set(gca, 'Ytick', y_values, 'YTickLabel',ylabels);
        title(sprintf('Delay Spread %3.0f ns',delay_spread(kk-skips)*1e9))
        ylabel(sprintf('Rx Pairs'))
        xlabel('TDOA (ns)')
        grid on
    end

    hleglines = [htrue(1) hcoarse(1) hrefine(1)];
    legend(ax(1), hleglines, 'True TDOA', 'Coarse TDOA Estimate', ...
        'Refined TDOA Estimate', 'Location', 'Northeast', ...
        'Position', [0.5 0.96 0.0 0.0])
    set(gcf, 'Position',  [100, 100, 1300, 800])
    
    annotation('textbox',[0.1 0.75 0.2 0.2],'String',my_anno, 'FitBoxToText','on');
end

%% Plot intermediate signals
if multi_option == 0 || multi_option == 4
    multi_choice = 1;
else
    multi_choice = 10;
end
trial_choice = 1;

figure
subplot(3,3,1)
plot(tx_raw{multi_choice}, '.-')
title('Tx Signal')
xlabel('Sample Number')
ylabel('Amplitude')
axis([-inf inf -inf inf])

subplot(3,3,2)
plot(tx_up{multi_choice}, '.-')
title('Tx Signal Upsampled')
xlabel('Sample Number')
ylabel('Amplitude')
axis([-inf inf -inf inf])

subplot(3,3,4)
base = 0;
for ii = 1:numrefs
    sig = tx_delayed{multi_choice}(:,ii);
    hc = max(abs(sig));
    if ii == 1
        plot(real(sig), '.-'); hold all
        base = base + 1.02*hc;
    else
        plot(real(sig) + base + hc, '.-');
        base = base + 2.02*hc;
    end
end
% for ii = 1:numrefs
%     plot(tx_delayed{multi_choice}(:,ii)+(ii-1)*4, '.-'); hold all
% end
title('Delay Added')
xlabel('Sample Number')
ylabel('Amplitude')
axis([-inf inf -inf inf])

% N = size(tx_delayed{multi_choice},1);
% tf = fft(tx_up{multi_choice}, N);
% rtf = fft(tx_delayed{multi_choice}, N);
% obj_fn1(targetPos, fc, fs, tx_pwr_dbm, refPos, rtf, tf)
% the actual attenuation due to path traveled is not implemented in this
% code. Think about how to fix the objective 

subplot(3,3,5)
base = 0;
for ii = 1:numrefs
    sig = tx_multi{multi_choice}{trial_choice}(:,ii);
    hc = max(abs(sig));
    if ii == 1
        plot(real(sig), '.-'); hold all
        plot(imag(sig), '.-')
        base = base + 1.02*hc;
    else
        plot(real(sig) + base + hc, '.-'); hold all
        plot(imag(sig) + base + hc, '.-')
        base = base + 2.02*hc;
    end
end
% for ii = 1:numrefs
%     plot(tx_multi{multi_choice}{trial_choice}(:,ii)+(ii-1)*4*max_nlos_amp, '.-'); hold all
% end
title('Multipath Added')
xlabel('Sample Number')
ylabel('Amplitude')
axis([-inf inf -inf inf])

subplot(3,3,7)
base = 0;
for ii = 1:numrefs
    sig = tx_noise{multi_choice}{trial_choice}(:,ii);
    hc = max(abs(sig));
    if ii == 1
        plot(real(sig), '.-'); hold all
        base = base + 1.02*hc;
        text(10,base,sprintf('%2.1f',avg_snr_db{multi_choice}(ii,trial_choice)))
    else
        plot(real(sig) + base + hc, '.-');
        base = base + 2.02*hc;
        text(10,base,sprintf('%2.1f',avg_snr_db{multi_choice}(ii,trial_choice)))
    end
end
title('Noise Added')
xlabel('Sample Number')
ylabel('Amplitude')
axis([-inf inf -inf 1.1*base])

% subplot(3,3,8)
% base = 0;
% for ii = 1:numrefs
%     sig = tx_agc{multi_choice}{trial_choice}(:,ii);
%     hc = max(abs(sig));
%     if ii == 1
%         plot(real(sig), '.-'); hold all
%         base = base + 1.02*hc;
%         text(10,base,sprintf('%2.1f',avg_snr_db{multi_choice}(ii,trial_choice)))
%     else
%         plot(real(sig) + base + hc, '.-');
%         base = base + 2.02*hc;
%         text(10,base,sprintf('%2.1f',avg_snr_db{multi_choice}(ii,trial_choice)))
%     end
%     
% end
% title('After AGC')
% xlabel('Sample Number')
% ylabel('Amplitude')
% axis([-inf inf -inf 1.1*base])

if ~contains(loc_alg_params.type, 'dpd')
    subplot(3,3,[3 6 9])
    base = 0;
    for ii = 1:numrefs-1
        sig = squeeze(corr_mag{multi_choice}(trial_choice,:,ii).');
        nsamps = length(sig);
        lags = (1:nsamps) - (nsamps+1)/2;
        hc = max(abs(sig));
        ylabels{ii} = sprintf('1,%i', ii+1);
        if ii == 1
            plot(lags/(fs*1e-9), real(sig), '.-'); hold all
            base = base + 1.0*hc;
            text(lags(1)/(fs*1e-9)+200, base,sprintf('%s',ylabels{ii}))
        else
            plot(lags/(fs*1e-9), real(sig) + base + hc, '.-');
            base = base + 2.0*hc;
            text(lags(1)/(fs*1e-9)+200, base,sprintf('%s',ylabels{ii}))
        end    
    end
    title('Correlator Outputs')
    xlabel('Lag (ns)')
    ylabel('Magnitude')
    axis([-inf inf -inf 1.1*base])
    grid on
end


%% Plot dpd objective function heatmap
if contains(loc_alg_params.type, 'dpd')% && ~contains(loc_alg_params.type, 'sparse')
    trial_choice = 1;
%     dpd_type = 'direct';
    dpd_type = 'vanilla';
    if strcmp(dpd_type, 'vanilla')
%         [xidx, yidx] = find(dpd_obj_vals{1,1}{trial_choice} == ...
%             max(max(dpd_obj_vals{1,1}{trial_choice})));
    else    
%         [xidx, yidx] = find(dpd_obj_vals{1,1}{trial_choice} == ...
%             min(dpd_obj_vals{1,1}{trial_choice}(:)));
        oldcmap = colormap;
        colormap( flipud(oldcmap) );
    end
    
    figure
    for kk = 1;%1:num_jj
        % plot objective function
        if strcmp(loc_alg_params.type,'sparse-dpd')
            hpc = pcolor(dpd_grid{1,1}{1}, dpd_grid{1,1}{2}, sparse_heatmap{1,1}{trial_choice}); hold all
        else
            hpc = pcolor(dpd_grid{1,1}{1}, dpd_grid{1,1}{2}, dpd_obj_vals{1,1}{trial_choice}); hold all
        end
        if strcmp(dpd_type, 'vanilla') || strcmp(loc_alg_params.type,'sparse-dpd')
            colormap('jet')
            h = colorbar;
            caxis('auto')
        else
            colormap('jet')
            oldcmap = colormap;
            colormap( flipud(oldcmap) );
            h = colorbar;
            mcolor = 1.2*median(dpd_obj_vals{1,1}{trial_choice}(:));
            caxis([0 mcolor]);
%             caxis('auto')
        end
        % shading interp;
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
        plot(targetPos(1),targetPos(2), 'gp', 'MarkerSize',10, 'MarkerFaceColor', 'g', 'LineWidth',1);
        
        % plot multipath reflection location
        if multi_option == 4
            for ii = 1:length(multi_coords)
                if ~isnan(multi_coords{ii})
                    for jj = 1:size(multi_coords{ii},2)
                        coord = multi_coords{ii}(:,jj);
                        plot(coord(1),coord(2), 'k^', 'MarkerSize',10, 'MarkerFaceColor', 'g');
                    end
                end
            end
            
            % add labels
            for ii = 1:length(multi_coords)
                if ~isnan(multi_coords{ii})
                    for jj = 1:size(multi_coords{ii},2)
                        coord = multi_coords{ii}(:,jj);
                        h = text(coord(1), coord(2), sprintf('%i', ii), ...
                            'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                            'FontSize', 8);
                        set(h, 'Color',[0, 0 ,0])
                    end
                end
            end
        end

        % plot dpd grid
        dpdh = plot(dpd_grid{1,1}{1}, dpd_grid{1,1}{2}, 'w.', 'MarkerFaceColor', 'k', ...
                'MarkerSize',4, 'HandleVisibility','off'); 
        legh(3) = dpdh(1);

%         % plot minimum obj val location
%         dpdh = plot(dpd_grid{1,1}{1}(xidx,yidx), dpd_grid{1,1}{2}(xidx,yidx), 'gd', ...
%                 'MarkerSize',8, 'linewidth', 1, 'MarkerFaceColor', 'k', ...
%                 'HandleVisibility','off'); 
            
        if strcmp(loc_alg_params.type,'sparse-dpd')
            % plot location of largest value of sparse vector
            dpdh = plot(sparse_max_coord{multi_choice}(1,trial_choice), ...
                sparse_max_coord{multi_choice}(2,trial_choice), 'go', ...
               'MarkerSize',5, 'linewidth', 1.3, 'MarkerFaceColor', 'k', ...
                    'HandleVisibility','off');
        else
            % plot minimum obj val location
            dpdh = plot(coords{multi_choice}(1,trial_choice), ...
                coords{multi_choice}(2,trial_choice), 'go', ...
               'MarkerSize',5, 'linewidth', 1.3, 'MarkerFaceColor', 'k', ...
                    'HandleVisibility','off');
        end
        
        if strcmp(loc_alg_params.type, 'ms-dpd1')
            % plot reflection location
            plot(rcoords{multi_choice}(1,trial_choice), rcoords{multi_choice}(2,trial_choice), 'wo', ...
               'MarkerSize',10, 'linewidth', 2, 'HandleVisibility','off');
        end
        
    end
end

%% Plot circles for each receiver los and nlos path
if multi_option == 4
    figure; hold on
    
    % Plot the LOS path range circles
    circ_colors = {'b','r','k','#77AC30','#4DBEEE','#0072BD','y'};
    for nn = 1:Nrx
        viscircles(refPos(:,nn).', ranges(nn), 'color', circ_colors{nn});
    end
    
    % Plot the NLOS path circles
    for ii = 1:length(multi_coords)
        if ~isnan(multi_coords{ii})
            for jj = 1:size(multi_coords{ii},2)
                mrange = c*multi_delays{ii}(:,jj);
                viscircles(refPos(:,ii).', mrange, 'color', ...
                    circ_colors{ii},'LineStyle','--');
            end
        end
    end
    axis tight
    
    if numrefs > 1
        % Plot the LOS path hyperbolas
        ref_only = 0;
        plot_multiple_hyperbolas(refPos, targetPos, ref_only, grid_bounds)
        axis square
        axis(grid_bounds)
    end
    
    % Plot the true toas and multipath toas
    xlims = xlim;
    ylims = ylim;
    a = 0;
    axis([xlims(1)-a xlims(2)+a ylims(1)-a ylims(2)+a+30])
    xlims = xlim;
    ylims = ylim;
    for ii = 1:numrefs
        mdelays = multi_delays{ii}*1e9;
        if isnan(mdelays)
            h = text(xlims(1)+2, 0.95*ylims(2)-0.05*(ii-1)*abs(ylims(2)), ...
            sprintf('%i: %3.1f', ii, toas_true(ii)*1e9), ...
            'horizontalalignment', 'left', 'verticalalignment', 'middle',...
            'FontSize', 8);
        else
            str_format = ['%i: %3.1f, ' repmat('%3.1f, ', 1, length(mdelays)-1) '%3.1f'];
            for jj = 1:length(mdelays)
                h = text(xlims(1)+2, 0.95*ylims(2)-0.05*(ii-1)*abs(ylims(2)), ...
                    sprintf(str_format, ii, toas_true(ii)*1e9, mdelays), ...
                    'horizontalalignment', 'left', 'verticalalignment', 'middle',...
                    'FontSize', 8);
            end
        end
        set(h, 'Color',[0, 0 ,0])
    end
    axis square
    
    % Plot the receiver locations
    numrefs = size(refPos,2);
    for ii = 1:numrefs
        plot(refPos(1,ii), refPos(2,ii), 'ks', 'MarkerFaceColor', 'k', ...
            'MarkerSize',10, 'HandleVisibility','off'); 
    end

    % Plot the receiver numbers
    for ii = 1:numrefs
        h = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle',...
            'FontSize', 8);
        set(h, 'Color',[1, 1 ,1])
    end
    
    % Plot the reflection locations
    for ii = 1:length(multi_coords)
        if ~isnan(multi_coords{ii})
            for jj = 1:size(multi_coords{ii},2)
                coord = multi_coords{ii}(:,jj);
                plot(coord(1),coord(2), 'k^', 'MarkerSize',10, 'MarkerFaceColor', 'g');
            end
        end
    end

    % Label the reflection locations
    for ii = 1:length(multi_coords)
        if ~isnan(multi_coords{ii})
            for jj = 1:size(multi_coords{ii},2)
                coord = multi_coords{ii}(:,jj);
                h = text(coord(1), coord(2), sprintf('%i', ii), ...
                    'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                    'FontSize', 8);
                set(h, 'Color',[0, 0 ,0])
            end
        end
    end
    
    % Plot the target location
    plot(targetPos(1),targetPos(2), 'kp', 'MarkerSize',10, 'MarkerFaceColor', 'g', 'LineWidth',1);
end

%% Plot hyperbolas
% if multi_option == 4 && numrefs > 1
if numrefs > 1
    figure; hold on
    
    % Plot the LOS path hyperbolas
    ref_only = 0;
    plot_multiple_hyperbolas(refPos, targetPos, ref_only, grid_bounds)
    axis square
    axis(grid_bounds)
    
    % Plot the NLOS hyperbolas
    
    % Plot the receiver locations
    numrefs = size(refPos,2);
    for ii = 1:numrefs
        plot(refPos(1,ii), refPos(2,ii), 'ks', 'MarkerFaceColor', 'k', ...
            'MarkerSize',10, 'HandleVisibility','off'); 
    end

    % Plot the true tdoas
    xlims = xlim;
    ylims = ylim;
    for ii = 1:numrefs-1
        h = text(xlims(1)+0.15*abs(xlims(1)), 0.95*ylims(2)-0.1*(ii-1)*abs(ylims(2)), ...
            sprintf('%i1: %3.1f', ii+1, tdoas_true(ii)*1e9), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle',...
            'FontSize', 8);
        set(h, 'Color',[0, 0 ,0])
    end
    
    % Plot the receiver numbers
    for ii = 1:numrefs
        h = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle',...
            'FontSize', 8);
        set(h, 'Color',[1, 1 ,1])
    end
    
    if multi_option == 4
        % Plot the reflection locations
        for ii = 1:length(multi_coords)
            if ~isnan(multi_coords{ii})
                for jj = 1:size(multi_coords{ii},2)
                    coord = multi_coords{ii}(:,jj);
                    plot(coord(1),coord(2), 'k^', 'MarkerSize',10, 'MarkerFaceColor', 'g');
                end
            end
        end

        % Label the reflection locations
        for ii = 1:length(multi_coords)
            if ~isnan(multi_coords{ii})
                for jj = 1:size(multi_coords{ii},2)
                    coord = multi_coords{ii}(:,jj);
                    h = text(coord(1), coord(2), sprintf('%i', ii), ...
                        'horizontalalignment', 'center', 'verticalalignment', 'middle',...
                        'FontSize', 8);
                    set(h, 'Color',[0, 0 ,0])
                end
            end
        end
    end
    
    % Plot the target location
    plot(targetPos(1),targetPos(2), 'kp', 'MarkerSize',10, 'MarkerFaceColor', 'g', 'LineWidth',1);
end

%% Plot the channel taps
% figure
% subplot(2,1,1)
% base = 0;
% for ii = 1:numrefs
%     sig = h{multi_choice}{trial_choice}(:,ii);
%     hc = max(abs(sig));
%     if ii == 1
%         plot(real(sig), '.-'); hold all
%         base = base + 1.02*hc;
%     else
%         plot(real(sig) + base + hc, '.-');
%         base = base + 2.02*hc;
%     end
%     
% end
% title('Real Channel Taps')
% xlabel('Tap Number')
% ylabel('Amplitude')
% axis([-inf inf -inf 1.2*base])
% 
% subplot(2,1,2)
% base = 0;
% for ii = 1:numrefs
%     sig = h{multi_choice}{trial_choice}(:,ii);
%     hc = max(abs(sig));
%     if ii == 1
%         plot(imag(sig), '.-'); hold all
%         base = base + 1.02*hc;
%     else
%         plot(imag(sig) + base + hc, '.-');
%         base = base + 2.02*hc;
%     end
%     
% end
% title('Imaginary Channel Taps')
% xlabel('Tap Number')
% ylabel('Amplitude')
% axis([-inf inf -inf 1.2*base])

%% Plot power delay profile
% figure
% num_rows = ceil(num_delay_spreads/4);
% for mm = 1:num_delay_spreads
%     subplot(num_rows,4,mm)
%     for nn = 1:Ntrials
%         stem(abs(h{1,mm}{nn})); hold on
%     end
%     title(sprintf('%3.0f ns', delay_spread(mm)*1e9))
% end


%%
% save('tdoa_test.mat')
% save('tdoa_data.mat', 'bx','by','stdx','stdy','coords','delay_spread', ...
%     'bias_coords','covar_coords','refPos','targetPos','num_jj','num_kk',...
%     'multi_option', 'mse_coords','num_delay_spreads', 'Tsym','tx_pwr_dbm',...
%     'corr_mag','fs','Ntrials','initial_coords','bounds','avg_snr_db')

% save('dpd_data.mat', 'bx','by','stdx','stdy','coords','delay_spread', ...
%     'bias_coords','covar_coords','refPos','targetPos','num_jj','num_kk',...
%     'multi_option', 'mse_coords','num_delay_spreads', 'Tsym','tx_pwr_dbm',...
%     'corr_mag','fs','Ntrials','initial_coords','bounds','avg_snr_db','grid_numx',...
%     'grid_numy','dpd_grid')
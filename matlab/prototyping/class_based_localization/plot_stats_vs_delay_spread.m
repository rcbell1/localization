clear; close all

% simulation parameters
sim_params.Ntrials = 100;   % number of independent noise/channel realizations

% receiver parameters
r1 = [22.8, -30*pi/180];
r2 = [22.8, -150*pi/180];
r3 = [22.8, 90*pi/180];
refPos = [ r1(1)*cos(r1(2)) r2(1)*cos(r2(2)) r3(1)*cos(r3(2)); ...   % equilateral triangle
           r1(1)*sin(r1(2)) r2(1)*sin(r2(2)) r3(1)*sin(r3(2))];
rx_params.locs = refPos;
rx_params.fs = 100e6;       % rx sample rate

% transmitter parameters
tx_params.type = 'prn';     % prn or ofdm
tx_params.fs = 50e6;        % tx sample rate
tx_params.loc = [sum(refPos(1,:))/3; sum(refPos(2,:))/3];
tx_params.pwr_dbm = 0;
tx_params.Nsym = 100;
tx_params.span = 10;
tx_params.sps = 2;
tx_params.symbol_rate = tx_params.fs/tx_params.sps;
tx_params.symbol_period = 1/tx_params.symbol_rate;
tx_params.excess_bw = 0.4;
tx_params.center_freq = 2.395e9;
            
% channel parameters
channel_params.delay_spread = 300e-9;   % time between first and last arrival
channel_params.multi_option = 0;    % 0,1,2 or 3
channel_params.num_paths = 2;       % option 1 only
channel_params.multi_idx_jump = 4;  % option 1 only
channel_params.num_delay_spreads = 10;  % option 2,3 only
channel_params.max_num_paths = inf;     % option 2 only

% localization algorithm params
loc_alg_params.type = 'dpd';

% instantiate the object
dpd = Localization(sim_params, tx_params, rx_params, channel_params, ...
    loc_alg_params);

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
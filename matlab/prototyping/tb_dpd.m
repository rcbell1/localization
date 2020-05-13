close all; clear

%% Needed by dpd
grid_xmin = -10;
grid_xmax = 10;
grid_ymin = -10;
grid_ymax = 10;
grid_numx = 100;
grid_numy = 100;

grid_def = [grid_xmin grid_xmax;
            grid_ymin grid_ymax;
            grid_numx grid_numy];
        
%% Receiver coords using polar coords
r1 = [22.8, -30*pi/180];
r2 = [22.8, -150*pi/180];
r3 = [22.8, 90*pi/180];
center = r1;
refPos = [ r1(1)*cos(r1(2)) r2(1)*cos(r2(2)) r3(1)*cos(r3(2)); ...   % equilateral triangle
           r1(1)*sin(r1(2)) r2(1)*sin(r2(2)) r3(1)*sin(r3(2))];
center = [center(1)*cos(center(2));center(1)*sin(center(2))];
refPos = -center + refPos; % origin centered equilateral triangle       

%% The figures will be bounded by this region
bcenter = [sum(refPos(1,:))/3; sum(refPos(1,:))/3; sum(refPos(2,:))/3; sum(refPos(2,:))/3].';
% bounds = bcenter + [-3 3 -3 3;
%                   -3 3 -3 3;
%                   -3 3 -3 3
%                   -3 3 -3 3];

num_emitters = size(targetPos,2);
% num_emitters = 10;
temp = ones(2,2*num_emitters);
temp(1,:) = -temp(1,:);
temp2 = reshape(temp,4,[]).';
a = sqrt(2*r1(1)^2+4*r1(1)*cos(120*pi/180));
bounds = bcenter + 1.0*a*reshape(temp,4,[]).';
              
%% Emitter pulse properties
tx_pwr_dbm = 10;         % emitter transmit power in dBm (USRP max is 10 dBm)
fs_tx = 200e6/70;
Nsym = 1000;              % number of symbols in signals
span = 10;              % total length of shaping filter in symbols
sps = 4;                % samples per symbol at the transmitter
fsym = fs_tx/sps;             % symbol rate of transmitter (signal bandwidth)
beta = 0.4;             % excess bandwidth of tx pulse shaping filter
fc = 2.395e9;             % center frequency of transmitter

%% Receiver properties
fs = 200e6/22;                % receiver sample rates (Hz)
wlen = 2*ceil(fs/fsym)+1; % moving maximum window length in samples, odd number
nstds = 9;                % number of standard deviations to declare peak
percent_of_peak = 0.8;    % get the number of samples needed on either side 
                          % of correlation peaks for the peak value to drop 
                          % by this percent for use in super resolution

% General simulation properties
c = 299792458;          % speed of light m/s
Ntrials = 100;            % keep at 1 for this plot, only for heat maps
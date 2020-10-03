classdef Localization
    %LOCALIZATION Estimates the coordinates of RF emitters
    %   This class enables the simulation of localization algorithms to
    %   to compare performance. The simulation parameters and 
    %   localization algorithm are user defined via input parameters.
    %
    %   For example this class can:
    %   1) generates a transmit waveform (PRN Sequence, OFDM, etc.)
    %   2) add proper delays between emitters and each receiver
    %   3) add multipath (number of paths and delay spread user defined)
    %   4) add noise based on distance traveled and free space path loss
    %   5) add carrier frequency offset
    %   6) add sample frequency offset
    %   7) compute the emitter location given the received complex baseband
    %      samples for each receiver
    %   
    %   Several localization algorithms are available to use:
    %   1) direct position determination (DPD)
    %   2) time difference of arrival (TDOA)
    %   3) angle of arrival (AOA)
    %   4) received signal strength (RSS)
    %
    %   If a localization algorithm does not exist, it
    %   is easy to define and add it to the class.
    
    properties
        % sim properites
        Ntrials (1,1) double {mustBeInteger,mustBePositive} = 100;
        
        % transmitter properties
        tx_signal_type (1,1) string {mustBeMember(tx_signal_type,{'prn','ofdm'})} = 'prn';
        tx_samp_rate (1,1) double {mustBePositive} = 1;
        tx_loc double = zeros(2,1);
        tx_pwr_dbm (1,1) double = 0;
        tx_Nsym (1,1) double {mustBeInteger} = 100;
        tx_span (1,1) double {mustBeInteger} = 2;
        tx_sps (1,1) double {mustBeInteger} = 2;
        tx_symbol_rate (1,1) double {mustBePositive} = 1;
        tx_symbol_period (1,1) double {mustBePositive} = 1;
        tx_excess_bw (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(tx_excess_bw,1)} = 0.4;
        tx_center_freq (1,1) double = 2.395e9;
        
        % receiver properties
        rx_samp_rate (1,1) double {mustBePositive} = 1;
        rx_samps
        rx_locs = zeros(2,3);
        
        % channel properties
        multi_option (1,1) double {mustBeNonnegative,mustBeLessThan(multi_option,4)} = 0;
        delay_spread (1,1) double {mustBePositive} = 300e-9;
        num_paths (1,1) double {mustBeNonnegative,mustBeInteger} = 0;
        num_delay_spreads (1,1) double {mustBeNonnegative,mustBeInteger} = 0;
        multi_idx_jump (1,1) double {mustBeNonnegative,mustBeInteger}
        max_num_paths (1,1) double {mustBeNonnegative}
        
        % localizaion algorithm properties
        loc_alg_type (1,1) string {mustBeMember(loc_alg_type,...
            {'dpd','tdoa-lsq','tdoa-si','tdoa-ts'})} = 'dpd';
        
        % output properties
    end
    
    methods
        function obj = Localization(sim_params,tx_params,rx_params,...
                channel_params,loc_alg_params)
            % sim properites
            obj.Ntrials = sim_params.Ntrials;
            
            % transmitter properties
            obj.tx_signal_type = tx_params.type;
            obj.tx_samp_rate = tx_params.fs;
            obj.tx_loc = tx_params.loc;
            obj.tx_pwr_dbm = tx_params.pwr_dbm;
            obj.tx_Nsym = tx_params.Nsym;
            obj.tx_span = tx_params.span;
            obj.tx_sps = tx_params.sps;
            obj.tx_symbol_rate = tx_params.symbol_rate;
            obj.tx_symbol_period = 1/obj.tx_symbol_rate;
            obj.tx_excess_bw = tx_params.excess_bw;
            obj.tx_center_freq = tx_params.center_freq;
            
            % receiver properties
            obj.rx_samp_rate = rx_params.fs;
            obj.rx_locs = rx_params.locs;
            
            % channel properties
            obj.delay_spread = channel_params.delay_spread;
            obj.multi_option = channel_params.multi_option;
            obj.num_paths = channel_params.num_paths;
            obj.num_delay_spreads = channel_params.num_delay_spreads;
            obj.multi_idx_jump = channel_params.multi_idx_jump;
            obj.max_num_paths = channel_params.max_num_paths;
            
            % localizaion algorithm properties
        end
        
        function obj = localize(obj)
            %LOCALIZE Estimate the emitter coordinates
            %   Run the full simulation and return emitter coordinates
            obj.coords = 1;
        end
        
        function tx_samps = generate_tx_signal(obj)
            
        end
    end
end


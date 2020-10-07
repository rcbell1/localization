function [out, channel_coeffs] = add_multipath(samples, s_t, fc, fs, ranges, delay_spread, ...
    num_paths, max_num_paths, multi_idx, multi_options, max_nlos_amp, ...
    min_num_taps, multi_dist_based, multi_delays, tx_pwr_dbm, show_plots)
% for each receiver, we convolve with an L path channel model that has
% magnitude between the direct path attenuation and zero and uniform phase

% Option 0: uses the delay spread to define the total taps and then fills
% in num_paths of them with non zero values
% Option 1:  
% mult_option = 1;

% min_num_taps = 50;
% max_nlos_amp = 1;

% For any delay_spread, there should be at least min_num_taps multipath taps to
% choose from regardless of the receiver sample rate. To achieve this the
% samples might need to be upsampled, then multipath applied and then
% downsampled back to the original rate

c = 299792458;      % speed of light m/s
lambda = c/fc;
[num_samps, numrefs] = size(samples);

channel_coeffs = ones(1,numrefs);
if multi_options == 0
    out = samples;
%     channel_coeffs = ones(1,numrefs);
elseif multi_options == 1
    if multi_idx == 0
        out = samples; % no multipath
    elseif num_taps < min_num_taps
        fs_new = ceil(min_num_taps/delay_spread);
        num_taps = min_num_taps;
        channel_coeffs = zeros(num_taps, numrefs);
        [p,q] = rat(fs_new/fs); % without this resample can hang
        new_samples = resample(samples,p,q);
        for ii = 1:numrefs
            a = [rand;rand];
            h_vals(ii) = max_nlos_amp*a(1).*exp(1j*2*pi*a(2));

            h = zeros(num_taps,1);
            h(1) = 1; % direct path uneffected here
            h(multi_idx) = h_vals(ii);
            temp_out(:,ii) = conv(new_samples(:,ii),h);
            channel_coeffs(:,ii) = h;
        end
        out = resample(temp_out, q, p);
        out = out(1:num_samps,:); % preserve original length
        
    else

        num_taps = ceil(delay_spread*fs); % number of samples per delay spread interval
        channel_coeffs = zeros(num_taps, numrefs);
        for ii = 1:numrefs
            a = [rand;rand];
            h_vals(ii) = max_nlos_amp*a(1).*exp(1j*2*pi*a(2));

            h = zeros(num_taps,1);
            h(1) = 1; % direct path uneffected here
            h(multi_idx) = h_vals(ii);
            out(:,ii) = conv(samples(:,ii),h);
            channel_coeffs(:,ii) = h;
        end
        
    end
elseif multi_options == 2
    num_taps = ceil(delay_spread*fs); % number of samples per delay spread interval
%     channel_coeffs = zeros(num_taps, numrefs);  
    if num_taps == 0
        out = samples;
    elseif num_taps < min_num_taps
        fs_new = ceil(min_num_taps/delay_spread);
        num_taps = min_num_taps;
        channel_coeffs = zeros(num_taps, numrefs);
        [p,q] = rat(fs_new/fs); % without this resample can hang
        new_samples = resample(samples,p,q);
        for ii = 1:numrefs
            num_paths = randi([0 min(max_num_paths, num_taps)]);
            multi_idx = sort(randperm(num_taps, max(0,num_paths-1))); % num_paths-1 bc one path is the direct path
    %         path_idxs = sort(randperm(num_taps-1,num_paths-1).'+1);
    %             range_multi = c/fs*path_idxs;
    %         h_vals = lambda./(4*pi*range_multi).* ...
    %             exp(1j*2*pi*rand(num_paths-1,1));

            h_vals = max_nlos_amp*rand(num_paths-1,1).* ...
                exp(1j*2*pi*rand(num_paths-1,1));

            h = zeros(num_taps,1);
            h(multi_idx) = h_vals(:);

            h(1) = 1; % direct path uneffected here

            temp_out(:,ii) = conv(new_samples(:,ii),h);
            channel_coeffs(:,ii) = h;
    %         h_vals = nan; % don't care in this scenario
        end
        out = resample(temp_out, q, p);
        out = out(1:num_samps,:); % preserve original length
        
    else
        
        channel_coeffs = zeros(num_taps, numrefs);
        for ii = 1:numrefs
            num_paths = randi([0 min(max_num_paths, num_taps)]);
            multi_idx = sort(randperm(num_taps, max(0,num_paths-1))); % num_paths-1 bc one path is the direct path
    %         path_idxs = sort(randperm(num_taps-1,num_paths-1).'+1);
    %             range_multi = c/fs*path_idxs;
    %         h_vals = lambda./(4*pi*range_multi).* ...
    %             exp(1j*2*pi*rand(num_paths-1,1));

            h_vals = max_nlos_amp*rand(num_paths-1,1).* ...
                exp(1j*2*pi*rand(num_paths-1,1));

            h = zeros(num_taps,1);
            h(multi_idx) = h_vals(:);

            h(1) = 1; % direct path uneffected here

            temp_out(:,ii) = conv(new_samples(:,ii),h);
            channel_coeffs(:,ii) = h;
        end
    end
elseif multi_options == 3
    num_taps = ceil(delay_spread*fs); % number of samples per delay spread interval
%     channel_coeffs = zeros(num_taps, numrefs);
    if num_taps == 0
        out = samples;
    elseif num_taps < min_num_taps
        fs_new = ceil(min_num_taps/delay_spread);
        num_taps = min_num_taps;
        channel_coeffs = zeros(num_taps, numrefs);
        [p,q] = rat(fs_new/fs); % without this resample can hang
        new_samples = resample(samples,p,q);
        for ii = 1:numrefs
            multi_idx_vals = 2:num_taps;
            
            multi_idx = sort(multi_idx_vals(randi([1,...
                length(multi_idx_vals)], 1,num_paths-1)));
            
            if multi_dist_based == 0
                % multipath power chosen randomly within limits
                h_vals = max_nlos_amp*rand(num_paths-1,1).* ...
                    exp(1j*2*pi*rand(num_paths-1,1));
             
                % debug stuff    
                %multi_idx = floor(num_taps/6*ii);
                %h_vals = max_nlos_amp;
            else
                % multipath power based on distance traveled
                range_multi = c/fs_new*multi_idx;
                h_vals = lambda./(4*pi*range_multi).* ...
                    exp(1j*2*pi*rand(num_paths-1,1));
            end
            


            h = zeros(num_taps,1);
            h(multi_idx) = h_vals(:);
            h(1) = 1; % direct path uneffected here

%             out(:,ii) = conv(samples(:,ii),h);
            temp_out(:,ii) = conv(new_samples(:,ii),h);
            channel_coeffs(:,ii) = h;
        end
        out = resample(temp_out, q, p);
        out = out(1:num_samps,:); % preserve original length
    else
        
        for ii = 1:numrefs
    %         num_paths = randi([0 min(max_num_paths, num_taps)]);
            multi_idx_vals = 2:num_taps;
%             multi_idx = sort(multi_idx_vals(randi([1,length(multi_idx_vals)], 1,num_paths-1)));
            
            multi_idx = sort(randperm(num_taps, max(0,num_paths-1))); % num_paths-1 bc one path is the direct path
            
            
            if multi_dist_based == 0
                % multipath power chosen randomly within limits
                h_vals = max_nlos_amp*rand(num_paths-1,1).* ...
                    exp(1j*2*pi*rand(num_paths-1,1));
             
                % debug stuff
                %multi_idx = floor(num_taps/6*ii);
                %h_vals = max_nlos_amp;
            else
                % multipath power based on distance traveled
                range_multi = c/fs*multi_idx;
                h_vals = lambda./(4*pi*range_multi).* ...
                    exp(1j*2*pi*rand(num_paths-1,1));        
            end

            h = zeros(num_taps,1);
            h(multi_idx) = h_vals(:);
            h(1) = 1; % direct path uneffected here

            out(:,ii) = conv(samples(:,ii),h);
            channel_coeffs(:,ii) = h;
%             h_vals = nan; % don't care in this scenario
        end
    end
    
elseif multi_options == 4
    % find the max delay to set the length of the output samples
%     tx_pwr = 10^(tx_pwr_dbm/10);
%     tx_amp = sqrt(tx_pwr);
    delays = [];
    for ii = 1:numrefs
        delays = [delays; multi_delays{ii}(:)];
    end
    max_delay = nanmax(delays);
    st_len = length(s_t);
    if ~isnan(max_delay)
        max_num_delay_samps = ceil(max_delay*fs);
        multi_len = max_num_delay_samps + st_len;
        out_len = max(multi_len, num_samps);
        out = zeros(out_len,numrefs);
        for ii = 1:numrefs
            if ~isnan(multi_delays{ii})
                delays = multi_delays{ii};
                curr_samps = samples(:,ii);
%                 curr_samps = s_t;
                if multi_dist_based == 0
                    mrval = max(abs(samples(:)));
                    msval = max(abs(s_t(:)));
                    mult = mrval/msval;
                    % random amplitude and phase
%                     multipaths = max_nlos_amp * rand(1,length(delays)) .* ...
%                         exp(1j*2*pi*rand(1,length(delays))) .* ...
%                         delayseq([samples(:,ii); ...
%                         zeros(max_num_delay_samps,1)], delays, fs);

                    % random phase
%                     multipaths = max_nlos_amp .* ...
%                         exp(1j*2*pi*rand(1,length(delays))) .* ...
%                         delayseq([samples(:,ii); ...
%                         zeros(max_num_delay_samps,1)], delays, fs);
                    
                    % amplitude and phase not random
%                     multipaths = max_nlos_amp .* ...
%                             delayseq([samples(:,ii); ...
%                             zeros(max_num_delay_samps,1)], delays, fs);
                        
                    multipaths = max_nlos_amp .* mult .* ...
                        delayseq([s_t; zeros(max_num_delay_samps,1)], ...
                        delays, fs);
                else
                    % with random phase
%                     multipaths = lambda./(4*pi*c*delays) .* ...
%                         exp(1j*2*pi*rand(1,length(delays))) .* ...
%                         delayseq([s_t; ...
%                         zeros(max_num_delay_samps,1)], delays, fs);
                    
                    % without random phase
                    multipaths = lambda./(4*pi*c*delays) .* ...
                        delayseq([s_t; zeros(max_num_delay_samps,1)], ...
                        delays, fs);
                end
                
%                 curr_samps(length(curr_samps)+max_num_delay_samps) = 0;
%                 out_len = max(length(multipaths), length(curr_samps));
                curr_samps(out_len) = 0;
                multipaths(out_len) = 0;
%                 out(out_len,:) = 0;
                out(:,ii)= curr_samps + sum(multipaths,2);
            else
                out(1:length(samples(:,ii)),ii) = samples(:,ii);
            end
        end
    else
        out = samples;
    end
else
    out = samples;
end

if show_plots == 1
    figure
    subplot(2,1,1)
    plot(real(samples)); hold on
    plot(imag(samples))
    title('Before Adding Multipath')
    xlabel('Samples')
    ylabel('Amplitude')
    subplot(2,1,2)
    plot(real(out)); hold on
    plot(imag(out))
    title('After Adding Multipath')
    xlabel('Samples')
    ylabel('Amplitude')
end
end


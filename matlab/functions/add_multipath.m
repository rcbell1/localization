function [out,h_vals] = add_multipath(samples, fc, fs, ranges, delay_spread, ...
    num_paths, max_num_paths, multi_idx, multi_options, show_plots)

% for each receiver, we convolve with an L path channel model that has
% magnitude between the direct path attenuation and zero and uniform phase

% Option 0: uses the delay spread to define the total taps and then fills
% in num_paths of them with non zero values
% Option 1:  
% mult_option = 1;

h_vals = nan;
if multi_options == 0
    out = samples;
elseif multi_options == 1
    if multi_idx == 0
        out = samples; % no multipath
        h_vals = nan;
    else
        [num_samps, numrefs] = size(samples);
        c = 299792458;      % speed of light m/s
%         lambda = c/fc;
        num_taps = ceil(delay_spread*fs); % number of samples per delay spread interval
        for ii = 1:numrefs
            a = [rand;rand];
            h_vals(ii) = 0.9*a(1).*exp(1j*2*pi*a(2));

            h = zeros(num_taps,1);
            h(1) = 1; % direct path uneffected here
            h(multi_idx+1) = h_vals(ii);
            out(:,ii) = conv(samples(:,ii),h);
        end

    end
elseif multi_options == 2
    [num_samps, numrefs] = size(samples);
    num_taps = ceil(delay_spread*fs); % number of samples per delay spread interval
            
    for ii = 1:numrefs
        num_paths = randi([0 min(max_num_paths, num_taps)]);
        multi_idx = sort(randperm(num_taps, max(0,num_paths-1))); % num_paths-1 bc one path is the direct path
%         path_idxs = sort(randperm(num_taps-1,num_paths-1).'+1);
%             range_multi = c/fs*path_idxs;
%         h_vals = lambda./(4*pi*range_multi).* ...
%             exp(1j*2*pi*rand(num_paths-1,1));
%         if num_paths-1 > 0
            h_vals = 0.9*rand(num_paths-1,1).* ...
                exp(1j*2*pi*rand(num_paths-1,1));

            h = zeros(num_taps,1);
            h(multi_idx) = h_vals(:);
%         end
        h(1) = 1; % direct path uneffected here
        
        out(:,ii) = conv(samples(:,ii),h);
        h_vals = nan; % don't care in this scenario
    end
elseif multi_options == 3
    [num_samps, numrefs] = size(samples);
    num_taps = ceil(delay_spread*fs); % number of samples per delay spread interval
    
    if num_taps == 0
        out = samples;
    else
%         pwr_multi_v_samps = 2;
%         max_avg_pwr = pwr_multi_v_samps*mean(abs(samples).^2);
        max_avg_pwr = 2;
        for ii = 1:numrefs
    %         num_paths = randi([0 min(max_num_paths, num_taps)]);
            multi_idx_vals = 2:num_taps;
            multi_idx = sort(multi_idx_vals(randi([1,length(multi_idx_vals)], 1,num_paths-1)));
%             multi_idx = sort(randperm(num_taps, max(0,num_paths-1))); % num_paths-1 bc one path is the direct path
            
            h_vals = sqrt(max_avg_pwr)*rand(num_paths-1,1).* ...
                exp(1j*2*pi*rand(num_paths-1,1));

            h = zeros(num_taps,1);
            h(multi_idx) = h_vals(:);
            h(1) = 1; % direct path uneffected here

            out(:,ii) = conv(samples(:,ii),h);
            h_vals = nan; % don't care in this scenario
        end
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


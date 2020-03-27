function [out] = add_multipath(samples, fc, fs, ranges, delay_spread, num_paths, show_plots)

% for each receiver, we convolve with an L path channel model that has
% magnitude between the direct path attenuation and zero and uniform phase
% the factor of 0.95 assures that any multipath is lower in received
% amplitude than the direct path as would be expected

if num_paths == 1
    out = samples;
else
    [num_samps, numrefs] = size(samples);
    c = 299792458;      % speed of light m/s
    lambda = c/fc;
    num_taps = ceil(delay_spread*fs); % number of samples per delay spread interval
    for ii = 1:numrefs
        path_idxs = sort(randperm(num_taps-1,num_paths-1).'+1);
        range_multi = c/fs*path_idxs;
%         h_vals = lambda./(4*pi*range_multi).* ...
%             exp(1j*2*pi*rand(num_paths-1,1));
        h_vals = rand(num_paths-1,1).* ...
            exp(1j*2*pi*rand(num_paths-1,1));

        h = zeros(num_taps,1);
        h(1) = 1; % direct path uneffected here
        h(path_idxs) = h_vals(:);
        out(:,ii) = conv(samples(:,ii),h);
    end
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


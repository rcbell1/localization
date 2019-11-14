function [out] = lower_samp_rate(x, down_rate, show_plots)

out = downsample(x, down_rate);

if show_plots == 1
    figure
    for ii = 1:size(x,2)
        plot(out(:,ii)+(ii-1)*1, '.-'); hold all
    end
    title('Downsampled To Receiver Rates')
    xlabel('Sample Number')
    ylabel('Amplitude')
end

end


function [out] = add_noise(x, tx_pwr_dbm, fc, ranges, show_plots)

c = 299792458;      % speed of light m/s
Gr = 0;             % receiver gain dBi, 0 for omnidirectional antenna
Gt = 0;             % transmitter gain dBi
thermal_pwr_dbm = -174; % thermal noise at room temp dBm
lambda = c/fc;      % wavelength of transmitted signal

[nrows, ncols] = size(x);

% free space path loss
rx_pwr_dbm = tx_pwr_dbm + Gr + Gt + 20*log10(lambda./(4*pi*ranges)); % dBm
snrdb = rx_pwr_dbm - thermal_pwr_dbm;
snr = 10.^(snrdb/10);

var_noise = 1./snr;        % noise power is variance
std_noise = sqrt(var_noise);

out = x + std_noise/sqrt(2).*(randn(nrows, ncols)+1j*randn(nrows,ncols));

if show_plots == 1
    figure
    for ii = 1:size(x,2)
        subplot(2,1,1)
        plot(real(out(:,ii)+(ii-1)*1), '.-'); hold all
        axis([-inf inf -1.5 6.5])
        title('Noise Added Real')
        xlabel('Sample Number')
        ylabel('Real Amplitude')
        subplot(2,1,2)
        plot(imag(out(:,ii)+(ii-1)*1j), '.-'); hold all
        axis([-inf inf -1 5])
        title('Noise Added Imaginary')
        xlabel('Sample Number')
        ylabel('Imaginary Amplitude')
    end

end

end


function [out] = add_noise(x, tx_pwr_dbm, noise_bw, fc, ranges, show_plots)

c = 299792458;      % speed of light m/s
Gr = 0;             % receiver gain dBi, 0 for omnidirectional antenna
Gt = 0;             % transmitter gain dBi
lambda = c/fc;      % wavelength of transmitted signal

[nrows, ncols] = size(x);

% noise power due to thermal noise and electronics
thermal_pwr_dbm = -174 + 10*log10(noise_bw); % thermal noise at room temp dBm
electronics_noise_figure_db = 10;   % account for electronics noise
noise_pwr_dbm = thermal_pwr_dbm + electronics_noise_figure_db;
noise_pwr = 10^(noise_pwr_dbm/10);
std_noise = sqrt(noise_pwr/2);

% free space path loss
rx_pwr_dbm = tx_pwr_dbm + Gr + Gt + 20*log10(lambda./(4*pi*ranges)); % dBm
rx_pwr = 10.^(rx_pwr_dbm/10);

x_norm = sqrt(rx_pwr).*x;

% snrdb = rx_pwr_dbm - noise_pwr_dbm;
% snr = 10.^(snrdb/10);

% var_noise = 1./snr;        % noise power is variance
% std_noise = sqrt(var_noise);
noise = std_noise.*(randn(nrows, ncols)+1j*randn(nrows,ncols));
out = x_norm + noise;

if show_plots == 1
    figure
    height = max(max(abs(out)));
    for ii = 1:size(x,2)
        subplot(2,1,1)
        plot(real(out(:,ii)+(ii-1)*height), '.-'); hold all
        axis([-inf inf -1.5 6.5])
        title('Noise Added Real')
        xlabel('Sample Number')
        ylabel('Real Amplitude')
        axis tight
        subplot(2,1,2)
        plot(imag(out(:,ii))+(ii-1)*height, '.-'); hold all
        axis([-inf inf -1 5])
        title('Noise Added Imaginary')
        xlabel('Sample Number')
        ylabel('Imaginary Amplitude')
        axis tight
    end
end

end


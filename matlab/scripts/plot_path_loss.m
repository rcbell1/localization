fc = 2.395e9;       % emitter center frequency
tx_pwr_dbm = 10;    % emitter transmit power dBm
Gr = 0;             % receiver gain dBi, 0 for omnidirectional antenna
Gt = 0;             % transmitter gain dBi

c = 299792458;      % speed of light m/s
lambda = c/fc;      % wavelength of transmitted signal
d = linspace(0,2*lambda,1000);   % distance from emitter to plot loss over

rx_pwr_dbm = tx_pwr_dbm + Gr + Gt + 20*log10(lambda./(4*pi*d));

figure
plot(d/lambda, rx_pwr_dbm, '.-'); hold all
plot([d(1) d(end)]/lambda, [tx_pwr_dbm tx_pwr_dbm], 'k--')
xlabel('d/\lambda')
ylabel('Received Power (dBm)')
legend('Rx Power','Tx Power')
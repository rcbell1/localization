function [out] = add_attenuation(x, fc, ranges)

c = 299792458;      % speed of light m/s
Gr = 0;             % receiver gain dBi, 0 for omnidirectional antenna
Gt = 0;             % transmitter gain dBi
lambda = c/fc;      % wavelength of transmitted signal

% [nrows, ncols] = size(x);

path_loss_dbm = Gr + Gt + 20*log10(lambda./(4*pi*ranges));
path_loss = 10.^(path_loss_dbm/10);
% path_loss = 1;

% If a test location lines up with a ref node range = 0 causes pwr to be
% inf. Change it to 1 representing no change to signal
inf_idx = isinf(path_loss);
if sum(inf_idx)
    path_loss(inf_idx) = 1;
end

out = sqrt(path_loss).*x;

end


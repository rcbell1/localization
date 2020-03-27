function [coords, grid, unique] = dpd(rx_samps, fs, refPos, grid_defs)
% This method uses the maximum eigenvalue over grid search method
% of Eq (15) in "Direct Position Determination of Narrowband Radio
% Frequency Transmitters" by Weiss
        
c = 299792458;
[N, M] = size(rx_samps);
N = 2^(nextpow2(N));

grid_xmin = grid_defs(1,1);
grid_xmax = grid_defs(1,2);
grid_ymin = grid_defs(2,1);
grid_ymax = grid_defs(2,2);
grid_numx = grid_defs(3,1);
grid_numy = grid_defs(3,2);
Ngrid = grid_numx*grid_numy;

gx = linspace(grid_xmin, grid_xmax, grid_numx);
gy = linspace(grid_ymin, grid_ymax, grid_numy);
[X,Y] = meshgrid(gx,gy);

% wk2 = 2*pi*(0:N-1).'*fs/N;
wk = 2*pi*ifftshift(((0:N-1)-N/2)).'*fs/N;
rx_f = fft(rx_samps,N);

p = zeros(2,Ngrid);
U = zeros(N,M);
max_evals = zeros(1,M);
for loc_idx = 1:Ngrid
    p(:,loc_idx) = [X(loc_idx);Y(loc_idx)]; % grid point coord
    delays = vecnorm(p(:,loc_idx)-refPos)/c; % time between p and each rx
    delays_diff = delays - delays(1);
%     a = % steering vector
    for t = 1:M
%         et2 = exp(1j*wk2.*delays_diff(t));
%         et = exp(1j*wk.*delays_diff(t));
%         U2(:,t) = exp(1j*wk2.*delays_diff(t)).*rx_f(:,t);
        U(:,t) = exp(1j*wk.*delays_diff(t)).*rx_f(:,t);
    end
    D = U'*U;
%     D2 = U2'*U2;
    max_evals(loc_idx) = max(eig(D));
%     max_evals2(loc_idx) = max(eig(D2));
end

[max_eval, max_eval_idx] = max(max_evals);
coords = p(:,max_eval_idx);
grid = p;
unique = 1;
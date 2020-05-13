function [coords, grid, objective_vals, unique] = dpd(rx_samps, fs, refPos, grid_defs)
% This method uses the maximum eigenvalue over grid search method
% of Eq (15) in "Direct Position Determination of Narrowband Radio
% Frequency Transmitters" by Weiss
        
c = 299792458;
[N, M] = size(rx_samps);
N = 2^(nextpow2(N));

% Create the search grid
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
[objX,objY] = meshgrid(1:grid_numx,1:grid_numy);

% wk2 = 2*pi*(0:N-1).'*fs/N;
wk = 2*pi*ifftshift(((0:N-1)-N/2)).'*fs/N;
rx_f = fft(rx_samps,N);

p = zeros(2,Ngrid);
U = zeros(N,M);
max_evals = zeros(1,M);
objective_vals = zeros(size(X));
for loc_idx = 1:Ngrid
    p(:,loc_idx) = [X(loc_idx);Y(loc_idx)]; % grid point coord
    delays = vecnorm(p(:,loc_idx)-refPos)/c; % time between p and each rx
    delays_diff = delays - delays(1);
%     a = % steering vector
    for t = 1:M
        U(:,t) = exp(1j*wk.*delays_diff(t)).*rx_f(:,t);
    end
    D = U'*U;
    max_evals(loc_idx) = max(eig(D));
    objective_vals(objY(loc_idx),objX(loc_idx)) = max_evals(loc_idx);
end

[max_eval, max_eval_idx] = max(max_evals);
objective_vals = objective_vals/max_eval;
coords = p(:,max_eval_idx);
% grid = p;
grid = {X,Y};
unique = 1;
end
function [coords, loc_grid, unique] = sparse_dpd(rx_t, fs, refPos, grid_defs)
% This version uses the toa approach but requires knowledge of the channel
% coefficients which is assumed known in this implementation

c = 299792458;
[Ns, Nr] = size(rx_t);
% N = 2 * 2^(nextpow2(N1));

% Create the emitter search grid
grid_xmin = grid_defs(1,1); grid_xmax = grid_defs(1,2);
grid_ymin = grid_defs(2,1); grid_ymax = grid_defs(2,2);
grid_numx = grid_defs(3,1); grid_numy = grid_defs(3,2);
Ng = grid_numx*grid_numy;
gx = linspace(grid_xmin, grid_xmax, grid_numx);
gy = linspace(grid_ymin, grid_ymax, grid_numy);

[X,Y] = meshgrid(gx,gy); % loop traverses vertically through grid
loc_grid = {X,Y};

wk = 2*pi*ifftshift(((0:Ns-1)-Ns/2)).'*fs/Ns;
% wk2 = 2*pi*(0:N-1).'*fs/N; % i don't understand why this is different

% create the A matrix for the sparse problem
rx_f = fft(rx_t,Ns);
A = zeros(Ns*Nr, Ng);
An = zeros(Ns*Nr, Ns);
r = rx_f(:); % stack received signal vectors
for loc_idx = 1:Ng % emitter grid index   

        p = [X(loc_idx);Y(loc_idx)]; % emitter location test coord    
        tau = vecnorm(p-refPos)/c; % ToF between p and each rx

        for ii = 1:Nr
            An(1+Ns*(ii-1):Ns*ii,:) = diag(exp(-1j*wk*tau(ii)));
        end
        sn = 1/Nr*An'*r;
        A(:,loc_idx) = An*sn;
end

% [ x_hat, support, runtime ] = L1min( A, r );
[ x_hat, support, runtime ] = L2rwmin( A, r );
[xidx,yidx] = ind2sub(size(X), support);

coords = [X(xidx,yidx); Y(xidx,yidx)];
unique = 1; % this does nothing but I need it for legacy reasons right now

end


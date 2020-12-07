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
gx = linspace(grid_xmin, grid_xmax, grid_numx)+1e-4;
gy = linspace(grid_ymin, grid_ymax, grid_numy)+1e-4;
[X,Y] = meshgrid(gx,gy);
[objX,objY] = meshgrid(1:grid_numx,1:grid_numy);

wk = 2*pi*ifftshift(((0:N-1)-N/2)).'*fs/N;
% wk2 = 2*pi*(0:N-1).'*fs/N;
% t1 = exp(1j*wk*1e-8);
% t2 = exp(1j*wk2*1e-8);
rx_f = fft(rx_samps,N);

p = zeros(2,Ngrid);
U = zeros(N,M);
max_evals = zeros(1,M);
objective_vals = zeros(size(X));
for loc_idx = 1:Ngrid
    p(:,loc_idx) = [X(loc_idx);Y(loc_idx)]; % grid point coord
    delays = vecnorm(p(:,loc_idx)-refPos)/c; % time between p and each rx
    delays_diff = delays - delays(1);
%     t3 = exp(1j*wk*delays_diff);
%     t4 = exp(1j*wk2*delays_diff);
%     a = % steering vector
    for t = 1:M
        U(:,t) = exp(1j*wk.*delays_diff(t)).*rx_f(:,t);
%         U2(:,t) = exp(1j*wk2.*delays_diff(t)).*rx_f(:,t);
    end
    D = U'*U;
    max_evals(loc_idx) = max(eig(D));
    objective_vals(objY(loc_idx),objX(loc_idx)) = max_evals(loc_idx);
    
    if loc_idx == coord_to_idx(13,10,grid_numx)
%         figure % time domains
%         for ii = 1:M
%             subplot(M,2,2*ii-1)
%             plot(real(rx_samps(:,ii)), 'b.-')
%             subplot(M,2,2*ii)
%             plot(imag(rx_samps(:,ii)), 'r.-')
%         end
%         figure % spectrums
%         for ii = 1:M
%             subplot(M,1,ii)
%             plot(20*log10(abs(rx_f(:,1))), 'b.-');
%         end
%             figure % spectrums
%             subplot(2,1,1)
%             plot(10*log10(tx_f), 'bx-'); hold on
%             plot(bin_idxs, 10*log10(tx_fm), 'ro')
%             subplot(2,1,2)
%             plot(20*log10(abs(rx_f(:,1))), 'bx-'); hold on
%             plot(20*log10(abs(rx_fm(:,1))), 'ro')
            
%             figure % estimated vs received signal
%             subplot(2,1,1)
%             fh(1) = plot(real(r), 'bx-'); hold on
%             fh(2) = plot(real(r_est), 'ro-');
%             title('Real Part Comparison')
%             xlabel('Sample Number')
%             ylabel('Amplitude')
%             ylims = ylim;
%             plot([Nsm Nsm],ylims, 'k--')
%             plot([2*Nsm 2*Nsm],ylims, 'k--')
%             legend(fh, 'Received', 'Estimated')
%             subplot(2,1,2)
%             fh(1) = plot(imag(r), 'bx-'); hold on
%             fh(2) = plot(imag(r_est), 'ro-');
%             title('Imaginary Part Comparison')
%             xlabel('Sample Number')
%             ylabel('Amplitude')
%             ylims = ylim;
%             plot([Nsm Nsm],ylims, 'k--')
%             plot([2*Nsm 2*Nsm],ylims, 'k--')
%             legend(fh, 'Received', 'Estimated')
            stop = 1;
        end
end

[max_eval, max_eval_idx] = max(max_evals);
objective_vals = objective_vals/max_eval;
coords = p(:,max_eval_idx);
% grid = p;
grid = {X,Y};
unique = 1;
end

function idx = coord_to_idx(row,col,Nrows)
    idx = (row-1)*Nrows + col;
end
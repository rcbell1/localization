function [coords, r_coord, loc_grid, reflect_grid, objective_vals, unique] = ...
    ms_dpd(rx_t, fs, fc, refPos, grid_defs, s_t, tx_pwr_dbm)
        

c = 299792458;
[N1, M] = size(rx_t);
Ndim = size(refPos,1);
N = 2 * 2^(nextpow2(N1));

wk = 2*pi*ifftshift(((0:N-1)-N/2)).'*fs/N;
% wk2 = 2*pi*(0:N-1).'*fs/N; % i don't understand why this is different
% t1 = exp(1j*wk*-1e-8);
% t2 = exp(1j*wk2*-1e-8);
lambda = c/fc;
tx_pwr = 10^(tx_pwr_dbm/10);
tx_amp = sqrt(tx_pwr);

% winr = kaiser(N1,6);
% winr = winr/sum(winr);
% wint = kaiser(length(s_t),6);
% wint = wint/sum(wint);
% rx_f = fft(rx_t.*winr,N);
% s_f = fft(s_t.*wint,N);
rx_f = fft(rx_t,N);
s_f = fft(s_t,N);

% Create the search grid
grid_xmin = grid_defs(1,1); grid_xmax = grid_defs(1,2);
grid_ymin = grid_defs(2,1); grid_ymax = grid_defs(2,2);
grid_numx = grid_defs(3,1); grid_numy = grid_defs(3,2);
Nlocgrid = grid_numx*grid_numy;
gx = linspace(grid_xmin, grid_xmax, grid_numx);
gy = linspace(grid_ymin, grid_ymax, grid_numy);

[X,Y] = meshgrid(gx,gy); % loop traverses vertically through grid
% [objX,objY] = meshgrid(1:grid_numx,1:grid_numy);
loc_grid = {X,Y};

% Create the reflection location grid
refgrid_numx = 5; refgrid_xmax= 0; refgrid_xmin = -10;
refgrid_numy = 5; refgrid_ymax = 35; refgrid_ymin = 25;
Nrefgrid = refgrid_numx*refgrid_numy;
gx = linspace(refgrid_xmin, refgrid_xmax, refgrid_numx);
gy = linspace(refgrid_ymin, refgrid_ymax, refgrid_numy);

[Xr,Yr] = meshgrid(gx,gy); % loop traverses vertically through grid
% [objX,objY] = meshgrid(1:grirefgrid_numxd_numx,1:refgrid_numy);
reflect_grid = {Xr,Yr};

obj_heatmap_vals = zeros(grid_numy,grid_numx);
objective_ref_vals = zeros(grid_numy,grid_numx,Nrefgrid+1);
for loc_idx = 1:Nlocgrid % emitter grid index

    p = [X(loc_idx);Y(loc_idx)]; % emitter location test coord

    Qij = 0;
    parfor ii = 2:M
        for jj = 1:(ii-1)
            p_i = refPos(:,ii); % ith receiver location
            p_j = refPos(:,jj); % jth receiver location
            r_i = rx_f(:,ii);   % ith receiver sample stream
            r_j = rx_f(:,jj);   % jth receiver sample stream

            r_ij = conj(r_i).*r_j;
            abs_s2 = abs(s_f).^2;

            Sij = zeros(N,3*Nrefgrid+1);
%             Sij = zeros(N,1);

            % The LOS only scenario
            dtau = (norm(p-p_i) - norm(p-p_j))/c;
            Sij(:,1) = abs_s2 .* exp(1j*wk*dtau);
%             Sij(:,1) = exp(1j*wk*dtau);
            
            if loc_idx == coord_to_idx(54,46,grid_numx)
                stop = 1;
            end

            % The multipath scenarios
            for ref_idx = 1:Nrefgrid % reflection scenario index
                rp = [Xr(ref_idx);Yr(ref_idx)]; % reflection test coord
                for kk = 1:3
                    switch kk
                        case 1
                            tau_i = norm(p-p_i);
                            tau_j = norm(p-rp)/c + norm(rp-p_j)/c;
                            dtau = tau_i - tau_j;
                        case 2
                            tau_i = norm(p-rp)/c + norm(rp-p_i)/c;
                            tau_j = norm(p-p_j);
                            dtau = tau_i - tau_j;
                        case 3
                            tau_i = norm(p-rp)/c + norm(rp-p_i)/c;
                            tau_j = norm(p-rp)/c + norm(rp-p_j)/c;
                            dtau = tau_i - tau_j;
                    end
%                         3*(ref_idx-1)+kk+1
                    Sij(:,3*(ref_idx-1)+kk+1) = abs_s2 .* exp(1j*wk*dtau);
%                     Sij(:,3*(ref_idx-1)+kk+1) = exp(1j*wk*dtau);
                    
                    if 3*(ref_idx-1)+kk+1 == 64
                        stop = 1;
                    end
                end
            end

            Qij = Qij + r_ij'*Sij*(Sij\r_ij);
        end
    end
    [row,col] = idx_to_coord(loc_idx,grid_numy);
    obj_heatmap_vals(row,col) = Qij;
end
objective_vals = real(obj_heatmap_vals/max(obj_heatmap_vals(:)));


% [mxv, midx] = min(objective_ref_vals(:));
% [~,~,p] = ind2sub(size(objective_ref_vals), midx);
% % [min_ref_val, min_ref_idx] = min(objective_ref_vals,[],3);
% objective_vals = objective_ref_vals(:,:,p);

% if p == 1
%     r_coord = nan;
% else
%     r_coord = [Xr(p-1);Yr(p-1)];
% end
r_coord = [nan;nan];

[xidx, yidx] = find(objective_vals == max(objective_vals(:)));

coords = [X(xidx,yidx); Y(xidx,yidx)];
unique = 1; % this does nothing but I need it for legacy reasons right now

end

% %             if loc_idx == 4794
% if loc_idx == coord_to_idx(54,46,grid_numx) && ...
%         ref_idx == coord_to_idx(3,3,refgrid_numx) %(3-1)*5 + 3
% %             if loc_idx == 3654 && ref_idx == 0
%     figure(1);
%     r_t2 = ifft(model_i);
%     r_t2 = r_t2(1:N1);
%     rd = abs(rx_t(:,ii)-r_t2);
%     mrd = max(rd);
%     subplot(1,M,ii)
%     plot(real(rx_t(:,ii)),'x-')
%     hold all
%     plot(real(r_t2),'o-')
% %                 plot(imag(r_t2),'--')
% %                 plot(rd, 'd-')
%     xlabel('Sample Number')
%     ylabel('Amplitude')
%     title(sprintf('Rx %i',ii))
%     legend('Rx Real','Grid Real')
% %                 legend('Rx Real','Grid Real','Grid Imag')
% end



function idx = coord_to_idx(row,col,Nrows)
    idx = (row-1)*Nrows + col;
end

function [row,col] = idx_to_coord(idx,Nx)
    row = mod(idx-1,Nx) + 1;
    col = ceil(idx/Nx);
end

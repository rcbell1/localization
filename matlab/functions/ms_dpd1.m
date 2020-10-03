function [coords, r_coord, loc_grid, reflect_grid, objective_vals, unique] = ...
    ms_dpd(rx_t, fs, fc, refPos, grid_defs, s_t, tx_pwr_dbm)
% This method uses the objective function given by Eq (4) in "Direct 
% Position Determination of Narrowband Radio Frequency Transmitters" by Weiss
% extended to include multipaths. To solve optimize the solution in the
% pressence of multipath a multipath reflection grid is used
        
% for test
% rx_t = rx_t(1:96,:);

c = 299792458;
[N1, M] = size(rx_t);
Ndim = size(refPos,1);
N = 2 * 2^(nextpow2(N1));
% rx_t(N,:) = 0; % increases frequency domain resolution to improve phase alignment

% This finds the actual solution to the objective function
% [coords, objective_vals1] = fmincon(@(p)obj_fn1(p, fc, fs, tx_pwr_dbm, refPos, rx_samps, tx_wf), p0);

% This is needed to plot the heatmap since fmincon only returns the final
% objective function value, not the value for every point on the grid we
% created. So this is like a manual grid search using the direct objective
% function
% function [obj_heatmap_vals, grid] = obj_fn_grid(grid_defs, fc, fs, tx_pwr_dbm, refPos, rx_t, s_t)
    
wk = 2*pi*ifftshift(((0:N-1)-N/2)).'*fs/N;
lambda = c/fc;
tx_pwr = 10^(tx_pwr_dbm/10);
tx_amp = sqrt(tx_pwr);

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
[objX,objY] = meshgrid(1:grid_numx,1:grid_numy);
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

model_order = 1; % how many reflection points to consider per receiver
obj_heatmap_vals = zeros(grid_numy,grid_numx);
objective_ref_vals = zeros(grid_numy,grid_numx,Nrefgrid+1);
for ref_idx = 0:Nrefgrid
    for loc_idx = 1:Nlocgrid

    %         p(:,loc_idx) = [X(loc_idx);Y(loc_idx)]; % grid point coord
        p = [X(loc_idx);Y(loc_idx)]; % emitter location test coord
        if ref_idx == 0
            rp = nan; % there will be no reflection
        else
            rp = [Xr(ref_idx);Yr(ref_idx)]; % reflection test coord
        end

        for ii = 1:M
            p_i = refPos(:,ii); % ith receiver location
            r_i = rx_f(:,ii);   % ith receiver sample stream
            tau_i = norm(p-p_i)/c; % LOS ToF
%             bi = tx_amp*lambda/(4*pi*c*tau_i); % LOS path loss
            bi = lambda/(4*pi*c*tau_i); % LOS path loss
%             los_model_i = s_f.*exp(-1j*wk*tau_i);
            los_model_i = bi.*s_f.*exp(-1j*wk*tau_i);
            
            if ref_idx == 0 || ii ~= 1
                tau_in = inf; % no reflection
                bin = 0;
                nlos_model_i = 0;
            else
                tau_in = norm(p-rp)/c + norm(rp-p_i)/c; % NLOS ToF
%                 bin = tx_amp*lambda/(4*pi*c*tau_in); % NLOS path loss
                bin = lambda/(4*pi*c*tau_in); % NLOS path loss
%                 nlos_model_i = 2*s_f.*exp(-1j*wk*tau_in);
                nlos_model_i = bin.*s_f.*exp(-1j*wk*tau_in);
            end          
            
%             nlos_model_i = 0;
            model_i = los_model_i + nlos_model_i;
            obj_val_i(ii) = norm(r_i - model_i)^2;
            obj_val_i_old(ii) = norm(r_i - los_model_i)^2;

%             if loc_idx == 4794
            if loc_idx == coord_to_idx(54,46,grid_numx) && ...
                    ref_idx == coord_to_idx(3,3,refgrid_numx) %(3-1)*5 + 3
%             if loc_idx == 3654 && ref_idx == 0
                figure(1);
                r_t2 = ifft(model_i);
                r_t2 = r_t2(1:N1);
                rd = abs(rx_t(:,ii)-r_t2);
                mrd = max(rd);
                subplot(1,M,ii)
                plot(real(rx_t(:,ii)),'x-')
                hold all
                plot(real(r_t2),'o-')
%                 plot(imag(r_t2),'--')
%                 plot(rd, 'd-')
                xlabel('Sample Number')
                ylabel('Amplitude')
                title(sprintf('Rx %i',ii))
                legend('Rx Real','Grid Real')
%                 legend('Rx Real','Grid Real','Grid Imag')
            end
        end
%         obj_heatmap_vals(objX(loc_idx),objY(loc_idx)) = sum(obj_val_i);
        [row,col] = idx_to_coord(loc_idx,grid_numy);
        obj_heatmap_vals(row,col) = sum(obj_val_i);
    end
%     objective_vals = obj_heatmap_vals/max(obj_heatmap_vals(:));
    objective_ref_vals(:,:,ref_idx+1) = obj_heatmap_vals;
end

[mxv, midx] = min(objective_ref_vals(:));
[~,~,p] = ind2sub(size(objective_ref_vals), midx);
% [min_ref_val, min_ref_idx] = min(objective_ref_vals,[],3);
objective_vals = objective_ref_vals(:,:,p);

if p == 1
    r_coord = nan;
else
    r_coord = [Xr(p-1);Yr(p-1)];
end

[xidx, yidx] = find(objective_vals == min(objective_vals(:)));

coords = [X(xidx,yidx); Y(xidx,yidx)];
unique = 1; % this does nothing but I need it for legacy reasons right now

end

function idx = coord_to_idx(row,col,Nrows)
    idx = (row-1)*Nrows + col;
end

function [row,col] = idx_to_coord(idx,Nx)
    row = mod(idx-1,Nx) + 1;
    col = ceil(idx/Nx);
end

function [coords, grid, objective_vals, unique] = ...
    dpd_direct(rx_samps, fs, fc, refPos, grid_defs, tx_wf, tx_pwr_dbm, p0)
% This method uses the objective function given by Eq (4) in "Direct 
% Position Determination of Narrowband Radio Frequency Transmitters" by Weiss
        
c = 299792458;
[N, M] = size(rx_samps);
N = 2^(nextpow2(N));
grid = nan;

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
% wk = 2*pi*ifftshift(((0:N-1)-N/2)).'*fs/N;
% rx_f = fft(rx_samps,N);

% p = zeros(2,Ngrid);
% U = zeros(N,M);
% max_evals = zeros(1,M);
% objective_vals = zeros(size(X));
% for loc_idx = 1:Ngrid
%     p(:,loc_idx) = [X(loc_idx);Y(loc_idx)]; % grid point coord
%     delays = vecnorm(p(:,loc_idx)-refPos)/c; % time between p and each rx
%     delays_diff = delays - delays(1);
% %     a = % steering vector
%     for t = 1:M
%         U(:,t) = exp(1j*wk.*delays_diff(t)).*rx_f(:,t);
%     end
%     D = U'*U;
%     max_evals(loc_idx) = max(eig(D));
%     objective_vals(objY(loc_idx),objX(loc_idx)) = max_evals(loc_idx);
% end
% 
% [max_eval, max_eval_idx] = max(max_evals);
% objective_vals = objective_vals/max_eval;
% coords = p(:,max_eval_idx);
% % grid = p;
% grid = {X,Y};


wk = 2*pi*ifftshift(((0:N-1)-N/2)).'*fs/N;
rx_f = fft(rx_samps,N);
s_f = fft(tx_wf,N);
lambda = c/fc;
% p0 = [-5;-5];
tx_pwr = 10^(tx_pwr_dbm/10);
tx_amp = sqrt(tx_pwr);

[coords, objective_vals1] = fmincon(@(p)obj_fn(p, lambda, tx_amp, wk, refPos, rx_f, s_f), p0);
[objective_vals, grid] = obj_fn_grid(grid_defs, lambda, tx_amp, wk, refPos, rx_f, s_f);
unique = 1;

end

function obj_val = obj_fn(p, lambda, tx_amp, wk, refPos, rx_f, s_f)

    M = size(refPos,2);
    c = 299792458;
    
    for ii = 1:M
        p_i = refPos(:,ii); 
        r_i = rx_f(:,ii);
        tau_i = norm(p-p_i)/c;
        bi = tx_amp*lambda/(4*pi*norm(p-p_i));
        obj_val_i(ii) = norm(r_i - bi.*s_f.*exp(-1j*wk*tau_i))^2;
%         fi = @(p) norm(rx_f - bi*s_f*exp(-1j*2*wk*tau_i))^2;
%         fi = @(p) norm(r_i - lambda/(4*pi*norm(p-p_i))*s_f*exp(-1j*2*wk*norm(p-p_i)/c))^2;
    end
    obj_val = sum(obj_val_i);
end

function [obj_heatmap_vals, grid] = obj_fn_grid(grid_defs, lambda, tx_amp, wk, refPos, rx_f, s_f)

    M = size(refPos,2);
    c = 299792458;
    
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
    grid = {X,Y};
    
    for loc_idx = 1:Ngrid
        p(:,loc_idx) = [X(loc_idx);Y(loc_idx)]; % grid point coord
        
        for ii = 1:M
            p_i = refPos(:,ii); 
            r_i = rx_f(:,ii);
            tau_i = norm(p-p_i)/c;
            bi = tx_amp*lambda/(4*pi*norm(p-p_i));
%             bi = abs(r_i);
            model_i = bi.*s_f.*exp(-1j*wk*tau_i);
            obj_val_i(ii) = norm(r_i - model_i)^2;
    %         fi = @(p) norm(rx_f - bi*s_f*exp(-1j*2*wk*tau_i))^2;
    %         fi = @(p) norm(r_i - lambda/(4*pi*norm(p-p_i))*s_f*exp(-1j*2*wk*norm(p-p_i)/c))^2;
        end
        obj_heatmap_vals(objY(loc_idx),objX(loc_idx)) = sum(obj_val_i);
%         obj_heatmap_vals(objX(loc_idx),objY(loc_idx)) = sum(obj_val_i);
    end
end
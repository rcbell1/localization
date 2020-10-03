function [coords, grid, objective_vals, unique] = ...
    dpd_direct(rx_samps, fs, fc, refPos, grid_defs, tx_wf, tx_pwr_dbm, p0)
% This method uses the objective function given by Eq (4) in "Direct 
% Position Determination of Narrowband Radio Frequency Transmitters" by Weiss
        
% c = 299792458;
% [N, M] = size(rx_samps);
% N = 2*2^(nextpow2(N));
% grid = nan;

% rx_f = fft(rx_samps,N);
% s_f = fft(tx_wf,N);

% Debug stuff
% rt = ifft(rx_f);
% st = ifft(s_f);
% st = st(1:length(tx_wf));
% drt = max(abs(rt-rx_samps));
% dst = max(abs(st-tx_wf));
% t=0;

% This finds the actual solution to the objective function
[coords, objective_vals1] = fmincon(@(p)obj_fn1(p, fc, fs, tx_pwr_dbm, refPos, rx_samps, tx_wf), p0);

% This is needed to plot the heatmap since fmincon only returns the final
% objective function value, not the value for every point on the grid we
% created. So this is like a manual grid search using the direct objective
% function
[objective_vals, grid] = obj_fn_grid(grid_defs, fc, fs, tx_pwr_dbm, refPos, rx_samps, tx_wf);

% Normalize the objective values to be between 0 and 1
% objective_vals = objective_vals/max(objective_vals(:));

unique = 1; % this does nothing but I need it for legacy reasons right now

end

% This version does not use RSS and assumes b_i = 1
function obj_val = obj_fn1(p, fc, fs, tx_pwr_dbm, refPos, rx_t, s_t)
    
    c = 299792458;
    [N,M] = size(rx_t);
    wk = 2*pi*ifftshift(((0:N-1)-N/2)).'*fs/N;
    lambda = c/fc;
    tx_pwr = 10^(tx_pwr_dbm/10);
    tx_amp = sqrt(tx_pwr);
    
    rx_f = fft(rx_t,N);
    s_f = fft(s_t,N);
    
%     figure
    for ii = 1:M
        p_i = refPos(:,ii); 
        r_i = rx_f(:,ii);
        tau_i = norm(p-p_i)/c;
%         bi = lambda/(4*pi*norm(p-p_i));
%         obj_val_i(ii) = norm(r_i - bi.*s_f.*exp(-1j*wk*tau_i))^2;
        obj_val_i(ii) = norm(r_i - s_f.*exp(-1j*wk*tau_i))^2;
        
%         r_t2 = ifft(s_f.*exp(-1j*wk*tau_i));
%         rd = abs(rx_t(:,ii)-r_t2);
%         mrd = max(rd);
%         subplot(1,M,ii)
%         plot(real(rx_t(:,ii)),'x-')
%         hold all
%         plot(real(r_t2),'.-')
%         plot(imag(r_t2),'o-')
%         plot(rd, 'd-')
    end
    obj_val = sum(obj_val_i);
end

% This version includes the RSS term
function obj_val = obj_fn2(p, fc, fs, tx_pwr_dbm, refPos, rx_t, s_t)

    c = 299792458;
    [N,M] = size(rx_t);
    wk = 2*pi*ifftshift(((0:N-1)-N/2)).'*fs/N;
    lambda = c/fc;
    tx_pwr = 10^(tx_pwr_dbm/10);
    tx_amp = sqrt(tx_pwr);
    
    rx_f = fft(rx_t,N);
    s_f = fft(s_t,N);
    
    for ii = 1:M
        p_i = refPos(:,ii); 
        r_i = rx_f(:,ii);
        tau_i = norm(p-p_i)/c;
        bi = lambda/(4*pi*norm(p-p_i));
        obj_val_i(ii) = norm(r_i - bi.*s_f.*exp(-1j*wk*tau_i))^2;
    end
    obj_val = sum(obj_val_i);
end

% This function is used to create the heatmap
function [obj_heatmap_vals, grid] = obj_fn_grid(grid_defs, fc, fs, tx_pwr_dbm, refPos, rx_t, s_t)
    
    c = 299792458;
    [N1,M] = size(rx_t);
    N = 2*2^(nextpow2(N1));
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
    Ngrid = grid_numx*grid_numy;
    gx = linspace(grid_xmin, grid_xmax, grid_numx);
    gy = linspace(grid_ymin, grid_ymax, grid_numy);
    
    [Y,X] = meshgrid(gx,gy); % loop traverses vertically through grid
    [objX,objY] = meshgrid(1:grid_numx,1:grid_numy);
    grid = {X,Y};
    
    for loc_idx = 1:Ngrid
%         p(:,loc_idx) = [X(loc_idx);Y(loc_idx)]; % grid point coord
        p = [X(loc_idx);Y(loc_idx)]; % test coord
        
        for ii = 1:M
            p_i = refPos(:,ii); 
            r_i = rx_f(:,ii);
            tau_i = norm(p-p_i)/c;
            bi = lambda/(4*pi*norm(p-p_i));
%             bi = abs(r_i);
            model_i = bi.*s_f.*exp(-1j*wk*tau_i);
%             model_i = s_f.*exp(-1j*wk*tau_i);
            obj_val_i(ii) = norm(r_i - model_i)^2;
            
%             if loc_idx == 4794
            if loc_idx == 3654
%             if loc_idx == 1654
            figure(1);
            r_t2 = ifft(model_i);
            r_t2 = r_t2(1:N1);
            rd = abs(rx_t(:,ii)-r_t2);
            mrd = max(rd);
            subplot(1,M,ii)
            plot(real(rx_t(:,ii)),'x-')
            hold all
            plot(real(r_t2),'o-')
%             plot(imag(r_t2),'--')
%             plot(rd, 'd-')
            xlabel('Sample Number')
            ylabel('Amplitude')
            title(sprintf('Rx %i',ii))
            legend('Rx Real','Grid Real')
%             legend('Rx Real','Grid Real','Grid Imag')
            end
        end
        % [22;-2] -> [X(46);Y(40)] Ngrid = (40-1)*80+46
        obj_heatmap_vals(objY(loc_idx),objX(loc_idx)) = sum(obj_val_i);
    end
    obj_heatmap_vals = obj_heatmap_vals/max(obj_heatmap_vals(:));
end
function obj_val = obj_fn1(p, fc, fs, tx_pwr_dbm, refPos, rx_f, s_f)

    [N,M] = size(rx_f);
    c = 299792458;
    wk = 2*pi*ifftshift(((0:N-1)-N/2)).'*fs/N;
    lambda = c/fc;
    tx_pwr = 10^(tx_pwr_dbm/10);
    tx_amp = sqrt(tx_pwr);
    
    for ii = 1:M
        p_i = refPos(:,ii); 
        r_i = rx_f(:,ii);
        tau_i = norm(p-p_i)/c;
%         bi = tx_amp*lambda/(4*pi*norm(p-p_i));
%         obj_val_i(ii) = norm(r_i - bi.*s_f.*exp(-1j*wk*tau_i))^2;
        obj_val_i(ii) = norm(r_i - s_f.*exp(-1j*wk*tau_i))^2;
    end
    obj_val = sum(obj_val_i);
end
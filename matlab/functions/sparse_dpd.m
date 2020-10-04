function [max_coord, coords, loc_grid, sparse_heatmap, unique] = sparse_dpd(rx_t, tx_t, fs, refPos, grid_defs)
% This version uses the tdoa approach 

c = 299792458;
[Ns, Nr] = size(rx_t);
% N = 2 * 2^(nextpow2(N1));

% Create the emitter search grid
grid_xmin = grid_defs(1,1); grid_xmax = grid_defs(1,2);
grid_ymin = grid_defs(2,1); grid_ymax = grid_defs(2,2);
grid_numx = grid_defs(3,1); grid_numy = grid_defs(3,2);
Ng = grid_numx*grid_numy;
gx = linspace(grid_xmin, grid_xmax, grid_numx)+1e-4;
gy = linspace(grid_ymin, grid_ymax, grid_numy)+1e-4;

[X,Y] = meshgrid(gx,gy); % loop traverses vertically through grid
loc_grid = {X,Y};

wk = 2*pi*ifftshift(((0:Ns-1)-Ns/2)).'*fs/Ns;
% wk2 = 2*pi*(0:N-1).'*fs/N; % i don't understand why this is different

% find frequency mask
rx_f = fft(rx_t,Ns);
rx_fa = abs(rx_f(:,1)).^2; % using first receiver to form mask
tx_f = abs(fft(tx_t,Ns)).^2;
fthresh = 0.01*max(rx_fa);
bin_idxs = find(rx_fa > fthresh);

% apply the mask to samples 
Nsm = length(bin_idxs);
rx_fm = rx_f(bin_idxs,:);
tx_fm = tx_f(bin_idxs);
wkm = wk(bin_idxs);

% create the A matrix for the sparse problem
rij_f = conj(rx_fm(:,2:end)).*rx_fm(:,1); % form the cross correlations
A = zeros(Nsm*(Nr-1), Ng);
An = zeros(Nsm*(Nr-1), Nsm);
r = rij_f(:); % stack received signal vectors
S = 1/(tx_fm'*tx_fm*(Nr-1))*(tx_fm*tx_fm');
parfor loc_idx = 1:Ng % emitter grid index   

        p = [X(loc_idx);Y(loc_idx)]; % emitter location test coord    
        tau = vecnorm(p-refPos)/c; % ToF between p and each rx
        dtau = tau(:,2:end)-tau(:,1);

        temp = cell(Nr-1,1);
        for ii = 1:(Nr-1)
%             O = diag(exp(1j*wk*dtau(ii)));
%             An(1+Ns*(ii-1):Ns*ii,:) = O;
%             O = diag(exp(1j*wkm*dtau(ii)));
%             An(1+Nsm*(ii-1):Nsm*ii,:) = O;
            temp{ii} = diag(exp(1j*wkm*dtau(ii)));
        end
        An = vertcat(temp{:});
        A(:,loc_idx) = An*S*An'*r;
        r_est = A(:,loc_idx);
        
%         if loc_idx == coord_to_idx(38,38,grid_numx)
%             figure % time domains
%             subplot(2,1,1)
%             plot(real(tx_t), 'b-'); hold on
%             plot(real(rx_t(:,1)), 'r.-')
%             subplot(2,1,2)
%             plot(imag(tx_t), 'b-'); hold on
%             plot(imag(rx_t(:,1)), 'r.-')
%             figure % spectrums
%             subplot(2,1,1)
%             plot(10*log10(tx_f), 'bx-'); hold on
%             plot(bin_idxs, 10*log10(tx_fm), 'ro')
%             subplot(2,1,2)
%             plot(20*log10(abs(rx_f(:,1))), 'bx-'); hold on
%             plot(bin_idxs, 20*log10(abs(rx_fm(:,1))), 'ro')
%             figure % estimated vs received signal
%             subplot(2,1,1)
%             plot(real(r), 'bx-'); hold on
%             plot(real(r_est), 'r.-')
%             subplot(2,1,2)
%             plot(imag(r), 'bx-'); hold on
%             plot(imag(r_est), 'r.-')
%             stop = 1;
%         end
end

% normalize the columns of A
% Anorms = sqrt(sum(abs((A.^2))));
% A = bsxfun(@rdivide,A,Anorms);

% [ x_hat, support, runtime ] = L1min_sdpd( A, r );
[ x_hat, x_no_thresh, support, runtime ] = L2rwmin_sdpd( A, r );
% [ x_hat, support, runtime ] = sbl( A, r );
[xidx,yidx] = ind2sub(size(X), support);
% [xidx,yidx] = ind2sub(size(X), [3]);
% [xidx,yidx] = ind2sub(size(X), [1,2]);

if length(support) == 0
    coords = [nan;nan];
else
    for ii = 1:length(xidx)
        coords(:,ii) = [X(xidx(ii),yidx(ii)); Y(xidx(ii),yidx(ii))];
    end
    [~, max_idx] = max(x_no_thresh);
    [max_xidx, max_yidx] = ind2sub(size(X), max_idx);
    max_coord = [X(max_xidx,max_yidx); Y(max_xidx,max_yidx)];
end
sparse_heatmap = reshape(normalize(abs(x_no_thresh),'range'), length(gx), []);
unique = 1; % this does nothing but I need it for legacy reasons right now

end

function idx = coord_to_idx(row,col,Nrows)
    idx = (row-1)*Nrows + col;
end

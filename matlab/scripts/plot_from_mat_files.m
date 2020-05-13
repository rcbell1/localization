clear; close all

load('dpd_data.mat');
tdoa = load('tdoa_data.mat');

if multi_option == 1
% temp = horzcat(bias_coords{:});
% bx = temp(1,:);
% by = temp(2,:);
% temp = horzcat(covar_coords{:});
% for jj = 1:num_multi_idxs
%     temp = diag(covar_coords{1,jj});
%     stdx(jj) = sqrt(temp(1));
%     stdy(jj) = sqrt(temp(2));
% end

% Plot multipath index vs statistics
figure
subplot(311)
plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), sqrt(mse_coords(1,:)), '.-'); hold all
ylims = ylim;
ylims(2) = min(60,ylims(2));
plot([Tsym:Tsym:delay_spread; Tsym:Tsym:delay_spread]/1e-9, [ylims(1) ylims(2)], 'k--');
title(sprintf('Two Paths with Increasing Delay Between Them, Tx Power %3.0f dBm', tx_pwr_dbm))
xlabel('Delay Between 1st and 2nd Path (ns)')
ylabel('RMSE (m)')
grid on
axis([-inf inf ylims(1) ylims(2)])

subplot(312)
plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), bx, '.-'); hold all
plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), by, '.-')
ylims = ylim;
ylims(1) = max(-60,ylims(1));
ylims(2) = min(60,ylims(2));
plot([Tsym:Tsym:delay_spread; Tsym:Tsym:delay_spread]/1e-9, [ylims(1) ylims(2)], 'k--');
xlabel('Delay Between 1st and 2nd Path (ns)')
ylabel('Bias (m)')
legend('b_x','b_y')
grid on
axis([-inf inf ylims(1) ylims(2)])

subplot(313)
plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), stdx, '.-'); hold all
plot((1:num_multi_idxs)*multi_jump/(fs*1e-9), stdy, '.-')
ylims = ylim;
ylims(2) = min(60,ylims(2));
plot([Tsym:Tsym:delay_spread; Tsym:Tsym:delay_spread]/1e-9, [ylims(1) ylims(2)], 'k--');
xlabel('Delay Between 1st and 2nd Path (ns)')
ylabel('Standard Deviation (m)')
legend('\sigma_x','\sigma_y')
grid on
axis([-inf inf ylims(1) ylims(2)])

plot_opts = {'b.-', 'r.-', 'g.-', 'm.-'};
num_opts = length(plot_opts);
kk = 1;
figure
% title('Correlation Peaks across Path Delays')
for ii = 1:num_multi_idxs
	subplot(num_multi_idxs,1,ii)
    nsamps = length(corr_mag_sqs{ii});
    lags = (1:nsamps) - (nsamps+1)/2;
    plot(lags/(fs*1e-9), corr_mag_sqs{ii}, plot_opts{kk}); hold on
    axis([-3*delay_spread 3*delay_spread -inf inf]/1e-9)
    
    if kk == num_opts
        kk = 1;
    else
        kk = kk + 1;
    end
end
elseif multi_option == 2 || multi_option == 3
    
% temp = horzcat(bias_coords{:});
% bx = temp(1,:);
% by = temp(2,:);
% temp = horzcat(covar_coords{:});
% for jj = 1:num_delay_spreads
%     temp = diag(covar_coords{1,jj});
%     stdx(jj) = sqrt(temp(1));
%     stdy(jj) = sqrt(temp(2));
% end

figure
subplot(3,2,1:2)
plot(delay_spread/1e-9, sqrt(mse_coords(1,:)), '.-'); hold all
plot(delay_spread/1e-9, sqrt(tdoa.mse_coords(1,:)), '.-');
ylims = ylim;
ylims(2) = min(60,ylims(2));
plot([Tsym:Tsym:max(delay_spread); Tsym:Tsym:max(delay_spread)]/1e-9, [ylims(1) ylims(2)], 'k--');
title(sprintf(['Two Paths, Tx Power %2.0f dBm, Avg Rx SNR Rx1:%3.1f dB, Rx2:%3.1f dB, Rx3:%3.1f dB'], ...
    tx_pwr_dbm, avg_snr_db{1}(1,1), avg_snr_db{1}(2,1), avg_snr_db{1}(3,1)))
xlabel('Delay Spread (ns)')
ylabel('Location RMSE (m)')
grid on
axis([-inf inf ylims(1) ylims(2)])
legend('dpd','tdoa')

subplot(323)
lh(1) = plot(delay_spread/1e-9, bx, 'b.-'); hold all
% lh(2) = plot(delay_spread/1e-9, by, 'r.-'); hold all
% plot(delay_spread/1e-9, by, '.-')
lh(2) = plot(delay_spread/1e-9, tdoa.bx, 'rx-');
% lh(4) = plot(delay_spread/1e-9, tdoa.by, 'rx-.')
ylims = ylim;
ylims(1) = max(-60,ylims(1));
ylims(2) = min(60,ylims(2));
temp = plot([Tsym:Tsym:max(delay_spread); Tsym:Tsym:max(delay_spread)]/1e-9, [ylims(1) ylims(2)], 'k--');
lh(3) = temp(1);
xlabel('Delay Spread (ns)')
ylabel('Bias x(m)')
grid on
axis([-inf inf ylims(1) ylims(2)])
legend(lh, 'dpd', 'tdoa','Symbol Boundaries')

subplot(324)
% lh(1) = plot(delay_spread/1e-9, bx, 'b.-'); hold all
lh(1) = plot(delay_spread/1e-9, by, 'b.-'); hold all
% plot(delay_spread/1e-9, by, '.-')
% lh(3) = plot(delay_spread/1e-9, tdoa.bx, 'bx-.')
lh(2) = plot(delay_spread/1e-9, tdoa.by, 'rx-');
ylims = ylim;
ylims(1) = max(-60,ylims(1));
ylims(2) = min(60,ylims(2));
temp = plot([Tsym:Tsym:max(delay_spread); Tsym:Tsym:max(delay_spread)]/1e-9, [ylims(1) ylims(2)], 'k--');
lh(3) = temp(1);
xlabel('Delay Spread (ns)')
ylabel('Bias y(m)')
grid on
axis([-inf inf ylims(1) ylims(2)])
legend(lh, 'dpd', 'tdoa','Symbol Boundaries')

subplot(325)
lh(1) = plot(delay_spread/1e-9, real(stdx), 'b.-'); hold all
lh(2) = plot(delay_spread/1e-9, real(tdoa.stdx), 'rx-');
ylims = ylim;
ylims(2) = min(60,ylims(2));
temp = plot([Tsym:Tsym:max(delay_spread); Tsym:Tsym:max(delay_spread)]/1e-9, [ylims(1) ylims(2)], 'k--');
lh(3) = temp(1);
xlabel('Delay Spread (ns)')
ylabel('Standard Deviation (m)')
grid on
axis([-inf inf ylims(1) ylims(2)])
legend(lh, 'dpd', 'tdoa','Symbol Boundaries')

subplot(326)
lh(1) = plot(delay_spread/1e-9, real(stdy), 'b.-'); hold all
lh(2) = plot(delay_spread/1e-9, real(tdoa.stdy), 'rx-');
ylims = ylim;
ylims(2) = min(60,ylims(2));
temp = plot([Tsym:Tsym:max(delay_spread); Tsym:Tsym:max(delay_spread)]/1e-9, [ylims(1) ylims(2)], 'k--');
lh(3) = temp(1);
xlabel('Delay Spread (ns)')
ylabel('Standard Deviation (m)')
grid on
axis([-inf inf ylims(1) ylims(2)])
legend(lh, 'dpd', 'tdoa','Symbol Boundaries')

% plot_opts = {'b.-', 'r.-', 'g.-', 'm.-'};
% num_opts = length(plot_opts);
% kk = 1;
% figure
% for ii = 1:num_delay_spreads
% 	subplot(num_delay_spreads,1,ii)
%     nsamps = length(corr_mag_sqs{ii});
%     lags = (1:nsamps) - (nsamps+1)/2;
%     plot(lags/(fs*1e-9), corr_mag_sqs{ii}, plot_opts{kk}); hold on
%     axis([-1*delay_spread(end) 1*delay_spread(end) -inf inf]/1e-9)
%     
%     if kk == num_opts
%         kk = 1;
%     else
%         kk = kk + 1;
%     end
% end
end

%% Plot localization results
% figure
% max_cols = 4;
% num_rows = ceil(max(num_jj,num_kk)/max_cols);  % 4 plots across max
% for kk = 1:max(num_jj,num_kk)
% %     plot_idx = mod(kk,max_cols);
%     subplot(num_rows,max_cols,kk)
%     numrefs = size(refPos,2);
%     for ii = 1:numrefs
%         plot(refPos(1,ii), refPos(2,ii), 'ks', 'MarkerFaceColor', 'k', ...
%             'MarkerSize',10, 'HandleVisibility','off'); 
%         hold all
%     end
%     plot(refPos(1,1), refPos(2,1), 'ks', 'MarkerFaceColor', 'k')
%     for ii = 1:numrefs
%         h = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
%             'horizontalalignment', 'center', 'verticalalignment', 'middle',...
%             'FontSize', 8);
%         set(h, 'Color',[1, 1 ,1])
%     end
%     plot(targetPos(1),targetPos(2), 'kx', 'MarkerSize',18);
% 
%     % plot dpd grid
%     plot(dpd_grid{1,1}{1}, dpd_grid{1,1}{2}, 'k.', 'MarkerFaceColor', 'k', ...
%             'MarkerSize',6, 'HandleVisibility','off'); 
%     
%     for ii = 1:Ntrials
%         plot(coords{1,kk}(1,ii),coords{1,kk}(2,ii), 'b.', 'MarkerSize',...
%             12, 'HandleVisibility','off');
%     end
%     plot(coords{1,1}(1),coords{1,1}(2), 'b.', 'MarkerSize',12)
%     plot(initial_coords(1), initial_coords(2), 'ko')
%     axis equal
%     axis(bounds(1,:))
%     xlabel('x (m)')
%     ylabel('y (m)')
%     xlen = bounds(1,2) - bounds(1,1);
%     ylen = bounds(1,4) - bounds(1,3);
%     text(bounds(1,2)-0.35*xlen,bounds(1,4)-0.15*ylen, ...
%     sprintf('b_x: %3.2f (m)\nb_y: %3.2f (m)\n\\sigma_x: %3.2e (m)\n\\sigma_y: %3.2e (m)', ...
%     bias_coords{kk}(1), bias_coords{kk}(2), sqrt(covar_coords{kk}(1,1)), ...
%     sqrt(covar_coords{kk}(2,2))), 'fontsize', 8)
% end
% legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'Position', [0.5 0.9 0.1 0.1])

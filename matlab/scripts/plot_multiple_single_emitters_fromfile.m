clear; close all
addpath('../functions')

plot_toa_countours = 0;
show_plots = 0;         % show plots for debugging
show_circles = 0;       % plot circles centered on emitter to visualize tdoa
show_hyperbolas = 0;    % plot hyperbolas to visualize intersection point
file_paths = {'../../data/3/tx_center/rx1/rx_pulses_sliced.mat';
    '../../data/3/tx_base/rx1/rx_pulses_sliced.mat';
    '../../data/3/tx_side/rx1/rx_pulses_sliced.mat';
    '../../data/3/tx_opp_side/rx1/rx_pulses_sliced.mat'};
% bounds = [-4 4 -4 4;
%           -4 4 -4 4;
%           -4 4 -4 4
%           -50 50 -50 50];
      
% file_paths = {'../../data/3/wired_center/rx9/rx_pulses_sliced.mat'};

% Equilateral Triangle
a = 3.9624;     % length of one side of desired equilateral triangle
b = sqrt(3)*a/2;
refPos = [ 0  -a          -a/2; ...   % equilateral triangle
           0  0      sqrt(3)*a/2];
% center = [sum(refPos(1,:))/3; sum(refPos(2,:))/3];  
center = [0;0];
refPos = -center + refPos; % origin centered equilateral triangle
% refPos = [[0;0] refPos];    % add a ref node at the center
targetPos3 = [-3*a/4;sqrt(3)/4*a];   % chair
targetPos2 = [-a/2;0];               % base
targetPos4 = [-a/4;sqrt(3)/4*a];     % opp chair
targetPos1 = [sum(refPos(1,:))/3; sum(refPos(2,:))/3];     % center

targetPos = [targetPos1 targetPos2 targetPos3 targetPos4];

bcenter = [sum(refPos(1,:))/3; sum(refPos(1,:))/3; sum(refPos(2,:))/3; sum(refPos(2,:))/3].';
bounds = bcenter + [-50 50 -50 50;
                  -3 3 -3 3;
                  -50 50 -50 50
                  -3 3 -3 3];
% bounds = bcenter + [-3 3 -3 3];
      
% Receiver properties
wlen = 20;                % moving maximum window length in samples
nstds = 9;                % number of standard deviations to declare peak
percent_of_peak = 0.9;    % get the number of samples needed on either side 
                          % of correlation peaks for the peak value to drop 
                          % by this percent for use in super resolution

% General simulation properties
c = 299792458;          % speed of light m/s

%% Estimate the location
jj = 1; % for future work
numfiles = size(file_paths,1);
for ii = 1:numfiles
    toas_true(ii,:) = vecnorm(refPos - targetPos(:,ii))/c*1e9;
    tdoas_true = toas_true(:,2:end) - toas_true(:,1:end-1)
    
    [coords{ii,jj}, bias_coords{ii,jj}, covar_coords{ii,jj}, ...
        mse_coords(ii,jj), tdoas_coarse{ii,jj}, tdoas_refined{ii,jj}, ...
        prob_correlation(ii,jj), prob_detection(ii,jj), Ntrials(ii,jj), ...
        unique(ii,jj)] = get_multiple_single_emitters_fromfile(...
        file_paths{ii}, targetPos, refPos, wlen, nstds, ...
        percent_of_peak, show_plots);
end

%% Plots
% Plot localization results
figure
for kk = 1:numfiles
    subplot(1,numfiles,kk)
    numrefs = size(refPos,2);
    for ii = 1:numrefs
        plot(refPos(1,ii), refPos(2,ii), 'ks', 'MarkerFaceColor', 'k', ...
            'MarkerSize',10, 'HandleVisibility','off'); 
        hold all
    end
    plot(refPos(1,1), refPos(2,1), 'ks', 'MarkerFaceColor', 'k')
    for ii = 1:numrefs
        h = text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle',...
            'FontSize', 8);
        set(h, 'Color',[1, 1 ,1])
    end
    plot(targetPos(1,kk),targetPos(2,kk), 'bx', 'MarkerSize',5);

    for nn = 1:length(coords{kk,1})
        plot(coords{kk,1}(1,nn),coords{kk,1}(2,nn), 'b.', 'MarkerSize',12, 'HandleVisibility','off');
    end
    plot(coords{kk,1}(1),coords{kk,1}(2), 'b.', 'MarkerSize',12)
    axis equal
    axis(bounds(kk,:))
    xlabel('x (m)')
    ylabel('y (m)')
    xlen = bounds(kk,2) - bounds(kk,1);
    ylen = bounds(kk,4) - bounds(kk,3);
    text(bounds(kk,2)-0.35*xlen,bounds(kk,4)-0.15*ylen, ...
     sprintf('b_x: %3.2f (m)\nb_y: %3.2f (m)\n\\sigma_x: %3.2e (m)\n\\sigma_y: %3.2e (m)', ...
     bias_coords{kk}(1), bias_coords{kk}(2), sqrt(covar_coords{kk}(1,1)), ...
     sqrt(covar_coords{kk}(2,2))), 'fontsize', 8)
    
     if plot_toa_countours == 1
        % Plot time contours
        tbounds = 2*bounds(kk,:);
        steps = .006*max(tbounds);
        [tx,ty] = meshgrid(tbounds(1):steps:tbounds(2),tbounds(3):steps:tbounds(4));
        tx2 = tx - targetPos(1,kk);
        ty2 = ty - targetPos(2,kk);
        tz = sqrt(tx2.^2+ty2.^2)/c*1e9;
        contour(tx,ty,tz, 4:15, '--', 'showtext', 'on', 'HandleVisibility', 'off')
     end
end
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'Position', [0.5 0.9 0.1 0.1])

%% Create table of TDOAs
for kk = 1:numfiles
    figure
    Coarse = tdoas_coarse{kk,1}*1e9;
    Fine = tdoas_refined{kk,1}*1e9;
    T1 = table(Coarse, Fine);
    uitable('ColumnWidth',{108 108 108 108}, 'Data',T1{:,:},'ColumnName',{'Coarse 12','Coarse 13','Fine 12','Fine 13'},...
    'RowName',T1.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 0.85]);

    True = tdoas_true(kk,:);
    T2 = table(True);
    uitable('ColumnWidth',{108 108}, 'Data',T2{:,:},'ColumnName',{'True 12','True 13'},...
    'RowName',T2.Properties.RowNames,'Units', 'Normalized', 'Position',[0.2, 0.85, 0.4, 0.15]);
    % annotation('textbox', [0.2, 0.9, 0.1, 0.1], 'String', sprintf('%3.2f\n',tdoas_true));
end


% Plot with hyperbolas
% subplot(1,2,2)
% for ii = 1:numrefs
%     plot(refPos(1,ii), refPos(2,ii), 'rx','HandleVisibility','off', 'MarkerSize',14); 
%     hold all
%     text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
% end
% plot(refPos(1,1), refPos(2,1), 'rx')
% plot(targetPos(1),targetPos(2), 'bo', 'MarkerSize',8);
% for ii = 1:size(coords,2)
%     plot(coords(1,ii),coords(2,ii), 'b.', 'MarkerSize',12, 'HandleVisibility','off');
% end
% plot(coords(1,1),coords(2,1), 'b.', 'MarkerSize',12)
% axis equal
% axis(bounds)
% xlabel('x (m)')
% ylabel('y (m)')
% title('Localization Results')
% legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'Position', [0.69 0.06 0.1 0.1])
% grid on
% 
% if show_circles == 1
%     viscircles(repmat(targetPos, 1, numrefs)',vecnorm(targetPos-refPos,2,1),...
%     'color', 'b', 'linestyle', '--', 'linewidth', 0.5);
% end
% 
% if show_hyperbolas == 1
%     plot_multiple_hyperbolas(refPos, targetPos, tdoas*c, bounds)
% end







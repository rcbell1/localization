clear; close all
addpath('../functions')

show_plots = 1;         % show plots for debugging
show_circles = 1;       % plot circles centered on emitter to visualize tdoa
show_hyperbolas = 1;    % plot hyperbolas to visualize intersection point
file_path = '../../data/tx_center/tx_center_record.mat';
bounds = [-5 5 -5 5];
targetPos = [0;0];      % known position
% bounds = [-150 150 -150 150];
% refPos = [-50 50 0; ... % triangle
%           -50 -50 50];

% Equilateral Triangle
a = 3.9624;     % length of one side of desired equilateral triangle
b = sqrt(3)*a/2;
refPos = [ 0  a          a/2; ...   % equilateral triangle
           0  0      sqrt(3)*a/2];
center = [sum(refPos(1,:))/3; sum(refPos(2,:))/3];       
refPos = -center + refPos; % origin centered equilateral triangle
% refPos = [[0;0] refPos];    % add a ref node at the center

% Receiver properties
wlen = 20;                % moving maximum window length in samples
nstds = 9;                % number of standard deviations to declare peak
percent_of_peak = 0.8;    % get the number of samples needed on either side 
                          % of correlation peaks for the peak value to drop 
                          % by this percent for use in super resolution

% General simulation properties
c = 299792458;          % speed of light m/s

%% Estimate the location
[coords, tdoas] = get_single_emitter_fromfile(file_path, refPos, wlen, ...
    nstds, percent_of_peak, show_plots);

%% Plots
% Plot localization results
figure
subplot(1,2,1)
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
plot(targetPos(1),targetPos(2), 'bo', 'MarkerSize',8);
for ii = 1:size(coords,2)
    plot(coords(1,ii),coords(2,ii), 'b.', 'MarkerSize',12, 'HandleVisibility','off');
end
plot(coords(1,1),coords(2,1), 'b.', 'MarkerSize',12)
axis equal
axis(bounds)
xlabel('x (m)')
ylabel('y (m)')
title('Localization Results')
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'Position', [0.25 0.06 0.1 0.1])
% grid on

% Plot time contours
tbounds = 2*bounds;
steps = round(.1*tbounds);
[tx,ty] = meshgrid(tbounds(1):steps:tbounds(2),tbounds(3):steps:tbounds(4));
tx2 = tx - targetPos(1);
ty2 = ty - targetPos(2);
tz = sqrt(tx2.^2+ty2.^2)/c*1e9;
contour(tx,ty,tz, '--', 'showtext', 'on', 'HandleVisibility', 'off')

subplot(1,2,2)
for ii = 1:numrefs
    plot(refPos(1,ii), refPos(2,ii), 'rx','HandleVisibility','off', 'MarkerSize',14); 
    hold all
    text(refPos(1,ii), refPos(2,ii), sprintf('%i', ii), 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
end
plot(refPos(1,1), refPos(2,1), 'rx')
plot(targetPos(1),targetPos(2), 'bo', 'MarkerSize',8);
for ii = 1:size(coords,2)
    plot(coords(1,ii),coords(2,ii), 'b.', 'MarkerSize',12, 'HandleVisibility','off');
end
plot(coords(1,1),coords(2,1), 'b.', 'MarkerSize',12)
axis equal
axis(bounds)
xlabel('x (m)')
ylabel('y (m)')
title('Localization Results')
legend('Receiver Locations', 'Target Emitter Location', 'Estimated Target Location', 'Position', [0.69 0.06 0.1 0.1])
grid on

if show_circles == 1
    viscircles(repmat(targetPos, 1, numrefs)',vecnorm(targetPos-refPos,2,1),...
    'color', 'b', 'linestyle', '--', 'linewidth', 0.5);
end

if show_hyperbolas == 1
    plot_multiple_hyperbolas(refPos, targetPos, tdoas*c, bounds)
end







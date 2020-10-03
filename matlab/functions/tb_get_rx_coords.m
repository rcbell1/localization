clear; close all

Nrx = 13;
radius = 7;
origin = 5;

[coords, ocoord] = get_rx_coords(Nrx, origin, radius);

% Plots
plot(coords(1,:), coords(2,:), '^','MarkerFaceColor','k'); hold on
fbound = 1.5*radius;
axis([ocoord(1)-fbound ocoord(1)+fbound ocoord(2)-fbound ocoord(2)+fbound])
axis equal

% make a reference circle
angs = 2*pi*(0:1/360:1-1/360);
cpnts = radius*exp(-1j*angs) + (ocoord(1)+1j*ocoord(2));
plot(real(cpnts),imag(cpnts),'--')
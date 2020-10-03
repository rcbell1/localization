function [rx_coords, center_coord] = get_rx_coords(Nrx, origin, radius_len)
% Nrx is the number of receiver coordinates you want returned
% radius_len is the length of the radius of the circle receivers lie on
% origin is a flag that determines which coordinate will be the (0,0)

point_angles = 2*pi*(0:1/Nrx:(1-1/Nrx));
if Nrx == 3
    points = radius_len*exp(1j*point_angles)*exp(1j*pi/2);
elseif Nrx == 4
    points = radius_len*exp(1j*point_angles)*exp(1j*pi/4);
else
	points = radius_len*exp(1j*point_angles);
end
rx_coords = [real(points); imag(points)];

% default origin is the center of the circle
origin_coord = [0;0];
if origin ~= 0
    origin_coord = rx_coords(:,origin);
    rx_coords = rx_coords - origin_coord;
end

center_coord = -origin_coord;
end


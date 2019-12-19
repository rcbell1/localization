function [coords, unique] = geo_sphere_int(refPos, tdoas)
% This function implements the spherical interpolation algorithm presented
% in "The Spherical Interpolcation Method of Source Localization" by Smith
% and Abel. It requires at least four ref stations for 2D localization

c = 299792458;              % speed of light m/s
unique = 1;                 % is solution unique?

% translate the spatial origin to the first sensor location
refPos_t = refPos - refPos(:,1);

S = refPos_t(:,2:end).';
d = c*tdoas';
delta = vecnorm(S,2,2).^2 - d.^2;

% Using backslash is better than pinv numerically
numerator = d.'*delta - d.'*S*(S\delta);
denom = d.'*d - d.'*S*(S\d);
Rs = 0.5*numerator/denom;

coords = 0.5*(S\(delta - 2*Rs*d)); % parenthesis really important here!
coords(isnan(coords)) = inf;   % random big number

% Translate the origin back to original location
coords = coords + refPos(:,1);
end


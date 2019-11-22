function [coords, unique] = geo_ls(refPos, tdoas)
% This function implements the spherical interpolation algorithm presented
% in "The Spherical Interpolcation Method of Source Localization" by Smith
% and Abel

c = 299792458;              % speed of light m/s
unique = 1;                 % is solution unique?

% translate the spatial origin to the first sensor location. Add this back
% to the final solution
refPos_t = refPos - refPos(:,1);

S = refPos_t(:,2:end).';
d = c*tdoas';
delta = vecnorm(S,2,2).^2 - d.^2;

% Using backslash is better than pinv numerically
% T = eye(N)-S*pinv(S);
% numerator = d.'*T*delta;
% denom = d.'*T*d;
numerator = d.'*delta - d.'*S*(S\delta);
denom = d.'*d - d.'*S*(S\d);
Rs = 0.5*numerator/denom;

% paren = delta - 2*Rs*d;
% coords1 = 0.5*S\paren;

coords = 0.5*S\(delta - 2*Rs*d);

% Translate the origin back
coords = coords + refPos(:,1);

 
end


function [coords] = geo_least_squares(refPos, tdoas)

numrefs = size(refPos,2);
c = 299792458;                      % speed of light m/s

d = c*tdoas';
x1 = refPos(1,1);
y1 = refPos(2,1);
for ii = 2:numrefs
    A(:,ii-1) = refPos(:,1) - refPos(:,ii);
    b(ii-1,:) = 0.5*(x1^2 - refPos(1,ii)^2 + y1^2 - ...
        refPos(2,ii)^2 + d(ii-1).^2);    
end
A = A.';

w = A\b;
g = A\d;

% Coefficients of quadratic equation
a1 = g(1)^2+g(2)^2-1;
b1 = 2*(w(1)*g(1)+w(2)*g(2)-x1*g(1)-y1*g(2));
c1 = x1^2+y1^2+w(1)^2+w(2)^2-2*x1*w(1)-2*y1*w(2);

% There will be two roots. First keep only real positive roots so it represents
% a valid distance. Then keep only points that are within the convex hull 
% of the reference stations. If no solutions are within the hull, then just
% pick the first one
r1 = roots([a1 b1 c1]);
idx = r1>=0;    % distances must be non-zero
r1 = r1(idx);

ridx = imag(r1) == 0;   % find indices of real values
r1p1 = r1(ridx);          % remove complex values

if isempty(r1p1)
    coords = w; % no solution, use zero
elseif length(r1) == 1
    coords = w + r1p1*g;  % unique solution, use it
else
    for ii = 1:2
        coords(:,ii) = w + r1p1(ii)*g;
    end
    k = convhull(refPos(1,:),refPos(2,:));
    for ii = 1:2
        in(ii) = inpolygon(coords(1,ii),coords(2,ii), refPos(1,k), refPos(2,k));
    end
    idx = in == 1; 
    if length(idx) > 1
        coords = coords(:,1);   % two solutions, use the first
    else
        coords = w + r1p1(idx)*g; % two solutions, pick the one inside the convex hull
    end
end

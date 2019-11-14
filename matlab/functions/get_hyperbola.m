function [out1,out2,keep_p] = get_hyperbola(rangeDiff, refLoc, target, bounds, both)
% think about removing the need for target location by using the sign of
% the tdoa/range_diff to decide which side of hyperbola to keep

s1 = refLoc(:,1);   % station 1
s2 = refLoc(:,2);   % station 2
diff = s1-s2;
distFoci = norm(diff)/2;
center = (s1+s2)/2;
theta = atan2d(diff(2),diff(1));
if theta == 180
    theta = 0;
end

a = rangeDiff/2;          % denominator coefficient of x for baseline hyperbola
b = sqrt(distFoci^2-a^2); % denominator coefficient of y for baseline hyperbola
rotate = theta;           % desired rotation angle of hyperbolas from baseline
translate = center;       % desired center of hyperbolas

xmin = bounds(1);
xmax = bounds(2);
ymin = bounds(3);
ymax = bounds(4);

x = floor(10*xmin):0.2:floor(10*xmax);% expand calculation bounds
y = @(x,a,b) sqrt(b^2*(x.^2/a^2-1));    % horizontal baseline hyperbolas

% Create the baseline hyperbolas
x1 = 0:0.2:floor(10*xmax);
y1 = real(y(x1,a,b));
idx1 = find(abs(real(y1))<eps);
y1(idx1(1:end-1)) = nan;
x1 = [x1 x1];
y1 = [y1 -y1];
x2 = floor(10*xmin):0.2:0;
y2 = real(y(x2,a,b));
idx2 = find(abs(real(y2))<eps);
y2(idx2(2:end)) = nan;
x2 = [x2 x2];
y2 = [y2 -y2];

% Rotate the baseline hyperbolas
R = [cos(pi/180*rotate) -sin(pi/180*rotate); ...
    sin(pi/180*rotate) cos(pi/180*rotate)];
Yb1 = [x1; y1];
Yr1 = R*Yb1;
xr1 = Yr1(1,:);
yr1 = Yr1(2,:);

Yb2 = [x2; y2];
Yr2 = R*Yb2;
xr2 = Yr2(1,:);
yr2 = Yr2(2,:);

% Translate the rotated hyperbolas
xrt1 = xr1 + translate(1);
xrt2 = xr2 + translate(1);
yrt1 = yr1 + translate(2);
yrt2 = yr2 + translate(2);

Hp = real([xrt1;yrt1]); % remove zero imaginary component from some samples
Hm = real([xrt2;yrt2]);

fp = [1;0];     % foci positive side
fn = [-1;0];    % foci negative side

fpnew = R*fp+translate; % track them through transforms
fnnew = R*fn+translate;

keep_p = norm(fpnew-target) < norm(fnnew-target);   % decide which hyperbola to keep

if both == 1
    out1= Hp;
    out2 = Hm;
else
    if keep_p == 1
        out1 = Hp;
        out2 = [nan;nan];
    else
        out1 = [nan;nan];
        out2 = Hm;
    end
end

end


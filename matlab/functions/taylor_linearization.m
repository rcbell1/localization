function [coords, unique] = taylor_linearization(refPos, tdoas, ...
    initial_coords, num_iters)

if isempty(initial_coords)
    initial_coords = [-20;15];
end
if isempty(num_iters)
    num_iters = 10;
end

num_dims = 2; % two dimensional coordinate estimation
[num_samps, num_refs] = size(refPos);
c = 299792458;              % speed of light m/s
unique = 1;                 % is solution unique?

f = @(xi,x1,t) norm(xi-t)-norm(x1-t);
dfx = @(xi,x1,t) (x1(1)-t(1))/norm(x1-t) - (xi(1)-t(1))/norm(xi-t);
dfy = @(xi,x1,t) (x1(2)-t(2))/norm(x1-t) - (xi(2)-t(2))/norm(xi-t);

J = zeros(num_refs-1, num_dims);
fvec = zeros(num_refs-1,1);
prev_coords = initial_coords;
for ii = 1:num_iters 
    % compute Jacobian
    for jj = 2:num_refs
        J(jj-1,:) = [dfx(refPos(:,jj), refPos(:,1), prev_coords) ...
                   dfy(refPos(:,jj), refPos(:,1), prev_coords)]; 
               
        fvec(jj-1) = f(refPos(:,jj), refPos(:,1), prev_coords);
    end
    
    % take a step towards solution
%     diff = c*tdoas.'-fvec;
%     coords1 = prev_coords + pinv(J)*diff;
%     coords = prev_coords + J\diff;
    d = c*tdoas.';
    coords = prev_coords + J\(d-fvec);
    prev_coords = coords;
end

if isnan(sum(coords))
    coords = initial_coords;
end

end

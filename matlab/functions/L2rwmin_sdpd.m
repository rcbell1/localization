function [ x_hat, x_no_thresh, support, runtime ] = L2rwmin_sdpd( A, b )

fstart = tic;

[m,n] = size(A);
w = ones(n,1);
x_hat = zeros(n,1);
x_prov = x_hat;
lambda = 0.1; % originally 0.01
epsilon = 1e-4;
Niter = 100; % originally 100

for kk = 1:Niter
    
    Q = diag(1./w);
    x_prov = Q*A'*((A*Q*A' + lambda*eye(m))\b);

    w = 1./(x_prov.^2 + epsilon);
    epsilon = epsilon*0.9;

end

% support = find(abs(x_prov) > 0.0001);
support = find(normalize(abs(x_prov),'range') > 0.8); % originally 0.0001
x_no_thresh = x_prov;
x_hat(support) = x_prov(support);
runtime = toc(fstart);

end
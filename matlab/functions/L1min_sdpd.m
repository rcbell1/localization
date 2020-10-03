function [ x_hat, support, runtime ] = L1min_sdpd( A, b )
% requires the cvx package to be installed

fstart = tic;

[m,n] = size(A);
x_hat = zeros(n,1);
x_prov = x_hat;

cvx_begin
    variable x_prov(n);
    minimize( norm(x_prov,1) );
    subject to
    A*x_prov == b;
cvx_end

% support = find(abs(x_prov) > 0.0001);
support = find(normalize(abs(x_prov),'range') > 0.8);

x_hat(support) = x_prov(support);
runtime = toc(fstart);

end


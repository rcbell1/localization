clear; close all
%% The symbol rate is identified by the number of samples you can fit into
% the width of the equivalent continuous time box pulse. In this case,
% whatever values num_ones takes on, this is how many samples fit into one
% box pulse, so the symbol rate is 1/num_ones

num_ones = 10;
num_zeros_L = 50;
num_zeros_R = 50;
x = [zeros(num_zeros_L,1); ones(num_ones,1); zeros(num_zeros_R,1)];

Nfft = 2^14;
faxis = -0.5:1/Nfft:0.5-1/Nfft;

X1 = 1/num_ones*fftshift(abs(fft(x,Nfft)).^2);
X2 = 10*log10(X1);

figure
subplot(3,1,1)
plot(x,'.-')

subplot(3,1,2)
plot(faxis, X1)
axis([-inf inf 0 num_ones+1])
grid on

subplot(3,1,3)
plot(faxis, X2)
axis([-inf inf -20 10*log10(num_ones)+1])
grid on

%% This section shows the same plots but for a series of boxes
N = 100;
sps = 6;
x = 2*randi([0 1], N, 1)-1; % 1's and -1's
% x = randi([0 1], N, 1);   % 0's and 1's
x = upsample(x,sps);
x = filter(ones(sps,1), 1, x);

Nfft = 2^14;
faxis = -0.5:1/Nfft:0.5-1/Nfft;

X1 = 1/N*fftshift(abs(fft(x,Nfft)).^2);
X1_max = max(X1);
X2 = 10*log10(X1);

figure
subplot(3,1,1)
plot(x,'.-')

subplot(3,1,2)
plot(faxis, X1)
axis([-inf inf 0 X1_max+1])
grid on

subplot(3,1,3)
plot(faxis, X2)
axis([-inf inf -20 10*log10(X1_max)+1])
grid on
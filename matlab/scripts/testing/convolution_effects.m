clear; close all

N = 60;     % length of first sequence
M = 3;      % length of second sequence
D = 20;      % integer delay
mu = 0;
sig = 0.01;
amp = 1;

gauss = @(x,mu,sig,amp)amp*exp(-( ( (x-mu).^2 ) / ( 2*sig.^2 ) )).';

% x = gauss(linspace(-4*sig, 4*sig, N), mu, sig, amp);
x = kaiser(N, 20);
y1 = [zeros(D,1); x];

h = [1 zeros(1,20) 1.7];
y2 = conv(y1,h);

figure
subplot(3,1,1)
plot(x)
xlabel('Sample Number')
ylabel('Amplitude')
title('Input')
axis([0 N+M+D-2 -inf inf])

subplot(3,1,2)
plot(y1)
xlabel('Sample Number')
ylabel('Amplitude')
title('Delayed Output')
axis([0 N+M+D-2 -inf inf])

subplot(3,1,3)
plot(y2)
xlabel('Sample Number')
ylabel('Amplitude')
title('Multipath Delayed Output')
axis([0 N+M+D-2 -inf inf])


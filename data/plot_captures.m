clear; close all

x1 = read_complex_binary('rx1.dat');
x2 = read_complex_binary('rx2.dat');
x3 = read_complex_binary('rx3.dat');

figure
subplot(3,3,1)
plot(real(x1)); hold all
plot(imag(x1))
subplot(3,3,2)
plot(real(x2)); hold all
plot(imag(x2))
subplot(3,3,3)
plot(real(x3)); hold all
plot(imag(x3))

bounds = [933000 936000];

y1 = x1(bounds(1):bounds(2));
y2 = x2(bounds(1):bounds(2));
y3 = x3(bounds(1):bounds(2));

subplot(3,3,4:6)
plot(real(y1)); hold all
plot(real(y2))
plot(real(y3))

[c12, lags12] = xcorr(y1,y2);
[c13, lags13] = xcorr(y1,y3);

subplot(3,3,7:9)
plot(lags12, abs(c12))
axis([lags12(1) lags12(end) -inf inf])
xlabel('Lag12')
ylabel('Correlation Magnitude')
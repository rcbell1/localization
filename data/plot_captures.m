clear; close all

nzeros = 1000;
npulse = 1000;
ns = nzeros+npulse;

x1 = read_complex_binary('rx1.dat');
x2 = read_complex_binary('rx2.dat');
x3 = read_complex_binary('rx3.dat');

len_x1 = length(x1);
ncenter = floor(length(x1)/2);
bounds = [ncenter-floor(ns/2) ncenter+floor(ns/2)]; % focus on single pulse

figure
subplot(4,3,1)
plot(real(x1)); hold all
plot(imag(x1))
title('Rx 192.168.10.2')
subplot(4,3,2)
plot(real(x2)); hold all
plot(imag(x2))
title('Rx 192.168.10.4')
subplot(4,3,3)
plot(real(x3)); hold all
plot(imag(x3))
title('Rx 192.168.10.5')

y1 = x1(bounds(1):bounds(2));
y2 = x2(bounds(1):bounds(2));
y3 = x3(bounds(1):bounds(2));
len_y1 = length(y1);

subplot(4,3,4:6)
plot(real(y1)); hold all
plot(real(y2))
plot(real(y3))
title('Zoomed in Correlation Region')
axis tight

[c12, lags12] = xcorr(y1,y2);
[c13, lags13] = xcorr(y1,y3);
cmag12 = abs(c12);
cmag13 = abs(c13);
[mval12, midx12] = max(cmag12);
[mval13, midx13] = max(cmag13);
clen = 2*len_y1-1;  % each xi stream should be equal length assumption
mlagidx12 = midx12-median(midx12);
mlagidx13 = midx13-median(midx13);

subplot(4,3,7:9)
plot(lags12, abs(c12))
axis([lags12(1) lags12(end) -inf inf])
xlabel('Lag12')
ylabel('Correlation Magnitude')
title('Cross-correlation Rx1 and Rx2')
text(mlagidx12+10,mval12/2,sprintf('Lag Index of Peak: %i\nPeak Value: %f', mlagidx12, mval12))

subplot(4,3,10:12)
plot(lags13, abs(c13))
axis([lags13(1) lags13(end) -inf inf])
xlabel('Lag12')
ylabel('Correlation Magnitude')
title('Cross-correlation Rx1 and Rx3')
text(mlagidx13+10,mval13/2,sprintf('Lag Index of Peak: %i\nPeak Value: %f', mlagidx13, mval13))
function [h] = channel_est_fft(x, r)
% Suppose you are given the channel input sequence x and the channel 
% output sequence r and you want to recover the channel coefficients 
% (or sequence) h, where r = conv(x,h). This function does that for you

tap_thresh = 0.01;  % magnitude cutoff for channel tap length estimate
Nout = 2*length(x);
fftr = fft(r, Nout);
fftx = fft(x, Nout);
h_f = fftr./fftx;
h = ifft(h_f);
num_taps_est = find(abs(h)>tap_thresh,1,'last');
h = h(1:num_taps_est);

end


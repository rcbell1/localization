function [out, noise_bw] = generate_signal(Nsym, fsym, sps, span, beta, show_plots)

noise_bw = fsym*(beta+1);

x1 = 2*randi([0 1], Nsym, 1)-1;
x2 = upsample(x1,sps);
rrc = rcosdesign(beta, span, sps);
rrc = rrc.'/max(rrc);

out = conv(rrc, x2);
out = out/sqrt(mean(abs(out).^2)); % normalize power to 1
% out = filter(rrc,1,x2);

if show_plots == 1
    figure
    subplot(4,1,1)
    stem(x1, 'filled')
    axis([-inf inf -2 2])
    title('Emitter Symbol Values')
    xlabel('Symbol Number')
    ylabel('Amplitude')
    
    subplot(4,1,2)
    plot(rrc)
    axis([-inf inf -0.3 1.1])
    title('Emitter Pulse Shape')
    xlabel('Sample Number')
    ylabel('Amplitude')
    text(length(rrc)/5,0.5,sprintf('Root Raised Cosine\nBeta = %2.1f\nSpan = %i Symbols\nSps = %i', beta, span, sps));
    
    subplot(4,1,3)
    plot(out)
    axis([-inf inf -2 2])
    title('Shaped Symbol Stream')
    xlabel('Sample Number')
    ylabel('Amplitude')
    
    subplot(4,1,4)
    fs_m = fsym*sps/1e6; % MHz
    flen = length(rrc);
    w=kaiser(flen,8)';
    w=w/sum(w);
    fsym_m = fsym/1e6;
    Nfft = 2*2^nextpow2(length(w));
    faxis = sps*fsym_m*(-0.5:1/Nfft:0.5-1/Nfft);
    spec = abs(1/Nfft*fft(rrc.*w,Nfft));
    spec = spec/max(spec);
    spec = fftshift(10*log10(spec));
    plot(faxis,spec)
    axis([-inf inf -90 5])
    title('Symbol Stream Spectrum')
    xlabel('Frequency (MHz)')
    ylabel('|FFT|^2 (dB)')
    grid on
    
%     axes('position',[0.6 0.16 0.2 0.08])
%     plot(faxis,spec); hold on
%     axis([-1.2*fsym_g 1.2*fsym_g -60 5])
%     title('Zoom Baseband')
%     xlabel('Frequency (GHz)')
%     ylabel('|FFT|^2 (dB)')
%     plot([0.002 0.002], [-60 5], 'k--')
%     plot([-0.002 -0.002], [-60 5], 'k--')
end
end


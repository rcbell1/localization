clear; close all

file_paths = {'../../data/14/tx_center/rfs9/1/rx_pulses_sliced.mat';
    '../../data/14/tx_center/rfs9/2/rx_pulses_sliced.mat';
    '../../data/14/tx_center/rfs9/3/rx_pulses_sliced.mat'};
nfiles = size(file_paths, 1);
pnum = 66;       % pulse number
Nfft = 128;    % number of fft points for plotting
titles = {'Wired','Wireless','Not Connected'};

ym = -inf;
ym_f = -inf;
for fnum = 1:nfiles
    load(file_paths{fnum})

    y = yblock(bounds(pnum):bounds(pnum+1)-1,:);    
    ylen = size(y,1);
    Nfft = 2^(nextpow2(4*ylen));
    temp = max(max(abs([real(y) imag(y)])));
    ym = max(ym,temp);

    faxis = fs_rx*(-0.5:1/Nfft:0.5-1/Nfft);
    win = kaiser(ylen,8);
    y_f = 10*log10(fftshift(1/sum(win.^2)*(abs(fft(win.*y,Nfft)))).^2 );
    temp_f = max(max(y_f));
    ym_f = max(ym_f,temp_f);

    y1_f = y_f(:,1);
    y2_f = y_f(:,2);
    y3_f = y_f(:,3);

    figure
    subplot(3,3,1)
    plot(real(y(:,1)))
    title(sprintf('%s Rx1 Real', titles{fnum}))
    ylabel('Amplitude')
    xlabel('Sample Number')
    axis([-inf inf -ym ym])
    
    subplot(3,3,2)
    plot(imag(y(:,1)))
    title(sprintf('%s Rx1 Imag', titles{fnum}))
    ylabel('Amplitude')
    xlabel('Sample Number')
    axis([-inf inf -ym ym])
    
    subplot(3,3,3)
    plot(faxis/1e6, y1_f)
    title(sprintf('%s Spectrum', titles{fnum}))
    ylabel('|S|^2')
    xlabel('Frequency (MHz)')
    axis([-inf inf -120 ym_f])

    subplot(3,3,4)
    plot(real(y(:,2)))
    title('Rx2 Real')
    ylabel('Amplitude')
    xlabel('Sample Number')
    axis([-inf inf -ym ym])
    
    subplot(3,3,5)
    plot(imag(y(:,2)))
    title('Rx2 Imag')
    ylabel('Amplitude')
    xlabel('Sample Number')
    axis([-inf inf -ym ym])
    
    subplot(3,3,6)
    plot(faxis/1e6, y2_f)
    title(sprintf('%s Spectrum', titles{fnum}))
    ylabel('|S|^2')
    xlabel('Frequency (MHz)')
    axis([-inf inf -120 ym_f])

    subplot(3,3,7)
    plot(real(y(:,3)))
    title('Rx3 Real')
    ylabel('Amplitude')
    xlabel('Sample Number')
    axis([-inf inf -ym ym])
    
    subplot(3,3,8)
    plot(imag(y(:,3)))
    title('Rx3 Imag')
    ylabel('Amplitude')
    xlabel('Sample Number')
    axis([-inf inf -ym ym])
    
    subplot(3,3,9)
    plot(faxis/1e6, y3_f)
    title(sprintf('%s Spectrum', titles{fnum}))
    ylabel('|S|^2')
    xlabel('Frequency (MHz)')
    axis([-inf inf -120 ym_f])
end

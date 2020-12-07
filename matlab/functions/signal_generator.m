function [out, noise_bw] = signal_generator(mod_type, Nsym, M, sps, fc, ...
    fs, beta, span)
% This function generates various modulation types for testing

% M = 2; % modulation order

x = randi([0 M-1], Nsym, 1); % random input data signal

switch mod_type
    case 'fsk'
        freq_sep = 3e6; % freq separation (Hz)
        noise_bw = fs/sps; % estimate of noise bandwidth
        out = fskmod(x, M, freq_sep, sps, fs);
        out = out/sqrt(mean(abs(out).^2)); % normalize power to 1
        
    case 'msk'
        noise_bw = fs/sps; 
        out = mskmod(x, sps);
        out = out/sqrt(mean(abs(out).^2)); % normalize power to 1
        
    case 'gfsk'
        mod = comm.CPMModulator(...
            'ModulationOrder', M, ...
            'FrequencyPulse', 'Gaussian', ...
            'BandwidthTimeProduct', 0.35, ...
            'ModulationIndex', 1, ...
            'SamplesPerSymbol', sps);
          meanM = mean(0:M-1);
        
        out = mod(2*(x-meanM));
        out = out/sqrt(mean(abs(out).^2)); % normalize power to 1
        noise_bw = fs/sps; 
        
    case 'cpfsk'
        mod = comm.CPFSKModulator(...
        'ModulationOrder', M, ...
        'ModulationIndex', 0.5, ...
        'SamplesPerSymbol', sps);
        meanM = mean(0:M-1);
        
        out = mod(2*(x-meanM));
        out = out/sqrt(mean(abs(out).^2)); % normalize power to 1
        noise_bw = fs/sps; 
        
    case 'psk'
        noise_bw = fs/sps; 
        phase_ini = pi/M;  % initial phase of constellation
        y1 = pskmod(x, M, phase_ini, 'gray');
        y2 = upsample(y1,sps);
        rrc = rcosdesign(beta, span, sps);
        rrc = rrc.'/max(rrc);
        out = conv(rrc, y2);
        out = out/sqrt(mean(abs(out).^2)); % normalize power to 1
        
    case 'pam'
        noise_bw = fs/sps; 
        phase_ini = pi/4;
        y1 = pammod(x, M, phase_ini, 'gray');
        y2 = upsample(y1,sps);
        rrc = rcosdesign(beta, span, sps);
        rrc = rrc.'/max(rrc);
        out = conv(rrc, y2);
        out = out/sqrt(mean(abs(out).^2)); % normalize power to 1
        
    case 'qam'
        noise_bw = fs/sps*(1+beta); 
        y1 = qammod(x, M, 'gray');
        y2 = upsample(y1,sps);
        rrc = rcosdesign(beta, span, sps);
        rrc = rrc.'/max(rrc);
        out = conv(rrc, y2);
        out = out/sqrt(mean(abs(out).^2)); % normalize power to 1
        
    case 'ofdm'
        Nfft  = 64;
        cplen = 16;
        nullIdx  = [1:6 33 64-4:64]';
        pilotIdx = [12 26 40 54]';

        numDataCarrs = Nfft-length(nullIdx)-length(pilotIdx);
        dataIn = complex(randn(numDataCarrs,Nsym),randn(numDataCarrs,Nsym));
        pilots = repmat(pskmod((0:3).',4),1,Nsym);

        x = randi([0 M-1], numDataCarrs, Nsym); % random input data signal
        qamSym = qammod(x,M,'UnitAveragePower',true);
        out = ofdmmod(qamSym,Nfft,cplen,nullIdx,pilotIdx,pilots);
        out = out/sqrt(mean(abs(out).^2)); % normalize power to 1
        noise_bw = fs*numDataCarrs/Nfft;
        
    case 'ssb'
        out = ssbmod(x,0,fs);
        out = out/sqrt(mean(abs(out).^2)); % normalize power to 1
        noise_bw = fs; 
        
    case 'sin'
        out = sin(2*pi*1/sps*(0:sps*Nsym-1)).';
        out = out/sqrt(mean(abs(out).^2)); % normalize power to 1
        noise_bw = fs; 
        
end

% x2 = 2*x-1;
% x2 = upsample(x2,sps);
% rrc = rcosdesign(beta, span, sps);
% rrc = rrc.'/max(rrc);
% 
% out2 = conv(rrc, x2);
% % out2 = [out2; zeros(500,1)];
% out2 = out2/sqrt(mean(abs(out2).^2)); % normalize power to 1
% 
% figure
% subplot(2,1,1)
% plot(real(out), 'b.-');  hold on
% % plot(real(out2), 'bo-'); 
% xlabel('Sample Number')
% ylabel('Amplitude')
% subplot(2,1,2)
% plot(imag(out), 'r.-'); hold on
% % plot(imag(out2), 'ro-')
% xlabel('Sample Number')
% ylabel('Amplitude')

end


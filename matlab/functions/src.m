function [out] = src(in, t_idx, delay, atten)

rt = in(:,t_idx);   % template signal
[Nsamps, Nrx] = size(in);

max_samps = Nsamps+max(delay);
out = zeros(max_samps,Nrx);
out(1:Nsamps,1) = rt;  % the template is not changed
for ii = 2:Nrx
    rtd = [zeros(delay(ii-1),1);rt];
    Nrtd = length(rtd);  
    
    ri = in(:,ii);  % current receiver
    ri(Nrtd) = 0;   % zero pad to make lengths equal
    
    ri_nlos = ri - atten(ii-1)*rtd;
%     ri_nlos = 1/atten(ii-1)*ri - rtd;
    
%     ri_hat = ri - ri_nlos;
    ri_hat = 1/atten(ii-1)*(ri - ri_nlos);
    ri_hat(max_samps) = 0;
    out(:,ii) = ri_hat;
    
%     figure
%     subplot(5,1,1)
%     plot(real(rt)); hold all
%     title('Template')
%     axis([0 300 -inf inf])
%     subplot(5,1,2)
%     plot(real(ri))
%     title('Received Signal i')
%     axis([0 300 -inf inf])
%     subplot(5,1,3)
%     plot(real(rtd))
%     title('Template Delayed')
%     axis([0 300 -inf inf])
%     subplot(5,1,4)
%     plot(real(ri_nlos)); hold all
%     title('Isolated NLOS')
%     axis([0 300 -inf inf])
%     subplot(5,1,5)
%     plot(real(ri_hat))
%     title('Recovered Template on Rx i')
%     axis([0 300 -inf inf])
end




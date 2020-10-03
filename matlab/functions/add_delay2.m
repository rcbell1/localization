function [out, toas, tdoas, ranges] = add_delay2(x, targetPos, refPos, fs, show_plots)

numrefs = size(refPos,2);             % number of reference stations
numsamps = length(x);               % number of signal samples
% numpairs = size(refPos,2)-1;        % number of unique tdoas
c = 299792458;                      % speed of light m/s
Ts = 1/fs;

for ii = 1:numrefs
    ranges(ii) = norm(targetPos - refPos(:,ii));
end
toas = ranges/c;

tdoas = nan;
if numrefs > 1
    for ii = 2:numrefs
        rangeDiffs(ii-1) = norm(targetPos - refPos(:,ii)) - norm(targetPos - refPos(:,1));
    end
    tdoas = rangeDiffs/c;
end

max_num_delay_samps = ceil((max(toas))/Ts);
out = zeros(max_num_delay_samps+numsamps,numrefs);
for ii = 1:numrefs
    xd = delayseq([x; zeros(max_num_delay_samps,1)], toas(ii), fs);
    out(:,ii) = [xd; zeros((max_num_delay_samps+numsamps)-length(xd),1)];
end
% Looks like I can replace the for loop above with this single line but it
% needs to be tested
% out2 = delayseq([x; zeros(max_num_delay_samps,1)], toas, fs);

% Stripping any leading zeros
% [~, first_nonzero_idxs] = max(out~=0);
% first_nonzero_idx = min(first_nonzero_idxs);
% out = out(first_nonzero_idx:end, :);

if show_plots == 1
    figure
    for ii = 1:numrefs
        plot(out(:,ii)+(ii-1)*1.5); hold all
    end
    title('Delay Added Corresponding to Emitter Location')
    xlabel('Sample Number')
    ylabel('Amplitude')
    axis([-inf inf -inf inf])
end
end


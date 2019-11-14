function [out, delays, tdoas, ranges] = add_delay(x, targetPos, refPos, fs, show_plots)

numrefs = size(refPos,2);           % number of reference stations
numpairs = size(refPos,2)-1;        % number of unique tdoas
c = 299792458;                      % speed of light m/s

for ii = 1:numrefs
    ranges(ii) = norm(targetPos - refPos(:,ii));
end

for ii = 2:numrefs
    rangeDiffs(ii-1) = norm(targetPos - refPos(:,ii)) - norm(targetPos - refPos(:,1));
end

toas = ranges/c;
Ts = 1/fs;
numSampDelay = round(toas/Ts);

out = zeros(max(numSampDelay)+length(x), numrefs);
for ii = 1:numrefs
    numsamps = numSampDelay(ii)+length(x);
    out(1:numsamps,ii) = [zeros(numSampDelay(ii),1); x];
end

delays = numSampDelay;
tdoas = rangeDiffs/c;

if show_plots == 1
    figure
    for ii = 1:numrefs
        plot(out(:,ii)+(ii-1)*1); hold all
    end
    title('Delay Added Corresponding to Emitter Location')
    xlabel('Sample Number')
    ylabel('Amplitude')
end
end


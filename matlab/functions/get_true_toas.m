function [toas, tdoas] = get_true_toas(rx_coords, tx_coord)
% returns the true toas given receiver and emitter coordinates

[Ndims, Nrx] = size(rx_coords);
c = 299792458;                      % speed of light m/s

ranges = norm(tx_coord - rx_coords(:,1));
for ii = 2:Nrx
    ranges(ii) = norm(tx_coord - rx_coords(:,ii));
    rangeDiffs(ii-1) = ranges(ii) - ranges(1);
end
toas = ranges/c;
tdoas = rangeDiffs/c;

end


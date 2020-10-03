function [] = plot_multiple_hyperbolas(refPos, target, ref_only, bounds)

both = 0;               % plot both sides of hyperbola = 1

c = 299792458;          % speed of light m/s
nrefs = size(refPos,2); % number of receive stations

if ref_only == 0
    % compute the TDOA between all pairs of reference stations
%     h = [];
    pair_indxs = nchoosek(1:nrefs,2);   % all possible pairs themselves
    npairs = size(pair_indxs,1);        % total number of possible pairs
    for kk = 1:npairs
        ii = pair_indxs(kk,1);
        jj = pair_indxs(kk,2);
        s1 = refPos(:,ii);
        s2 = refPos(:,jj);
        tdoa = abs(norm(s1-target)-norm(s2-target))/c;
        [Hp,Hm] = get_hyperbola(tdoa*c, [s1 s2], target, bounds, both);
        plot_hyperbola(Hp, Hm);
    end
else
    % compute the TDOA between the reference and other receive stations
%     h = [];
    pair_indxs = nchoosek(1:nrefs,2);   % the unique pairs themselves
    npairs = nrefs - 1;                 % total number of unique pairs
    s1 = refPos(:,1);                   % the reference receiver
    for kk = 1:npairs
%         ii = pair_indxs(kk,1);
%         jj = pair_indxs(kk,2);
%         s1 = refPos(:,ii);
        s2 = refPos(:,kk+1);
        tdoa = abs(norm(s1-target)-norm(s2-target))/c;
        [Hp,Hm] = get_hyperbola(tdoa*c, [s1 s2], target, bounds, both);
        plot_hyperbola(Hp, Hm);
    end
end

end
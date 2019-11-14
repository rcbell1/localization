function [] = plot_hyperbola(Hp, Hm)

numHyp = size(Hp,1)/2;

for ii = 1:2:numHyp
    plot(Hp(ii,:), Hp(ii+1,:), 'HandleVisibility','off'); hold on
    plot(Hm(ii,:), Hm(ii+1,:), 'HandleVisibility','off')
end

% plot(Hp(1,:), Hp(2,:)); hold on
% plot(Hm(1,:), Hm(2,:))

end


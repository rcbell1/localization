function [out, fs, bounds] = load_signals_fromfile(file_path, show_plots)

load(file_path)

fs = fs_rx;
out = yblock;

if show_plots == 1
    up = max(max(real(out)));
    low = min(min(real(out)));
    delta = max(abs(up), abs(low));
    delta = 1.1*delta;
    for nn = 1:length(bounds)-1
        y = out(bounds(nn):bounds(nn+1)-1,:);
        plot(real(y) + delta*(nn-1)); hold on
%         plot(imag(y) + delta*(nn-1));
    end
end

end


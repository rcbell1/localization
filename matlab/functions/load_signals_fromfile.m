function [out, fs, bounds] = load_signals_fromfile(file_path, show_plots)

load(file_path)

fs = fs_rx;
out = yblock;

end


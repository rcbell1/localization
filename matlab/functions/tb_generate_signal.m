% testbench for generate_signal

Nsym = 50;  % total number of symbols
sps = 20;   % total samples per symbol
span = 20;  % total symbols the shaping filter spans

x = generate_signal(Nsym,sps,span);

plot(x)
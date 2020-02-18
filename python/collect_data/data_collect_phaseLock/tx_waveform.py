import numpy as np
import commpy as cp
import scipy.linalg as la
import matplotlib.pyplot as plt

# Construct TX waveform
def transmit_waveform(Nsym=100,sps=2,span=10,rolloff=0.5,pad_zeros=1500,plot_debug=False):
  # symbols are just the raw symbols with no shaping
  # tx_stream is the raw convolution of symbols with pulse reshape
  # tx_stream_clean removes the transient filter effects from the front and back of tx_stream

    num_samples_rrc = span*sps+1
    Ts = 1.0
    Fs = sps
    rrc = cp.filters.rrcosfilter(num_samples_rrc,rolloff,Ts,Fs)[1]
    rrc = rrc/(np.max(rrc))
    symbols = 2*np.random.randint(2,size=(Nsym))-1
    symbols_up = np.zeros(sps*Nsym)
    symbols_up[::sps] = symbols

    tx_stream = np.convolve(rrc,symbols_up)
    tx_stream_clean = tx_stream[sps*span/2:-(sps*span/2+1)]
    stream_len = tx_stream.size
    if(plot_debug):
      plt.plot(tx_stream)
      plt.plot(range(sps*span/2,tx_stream.size-sps*span/2-1), tx_stream_clean)
      plt.plot(range(sps*span/2,tx_stream.size-sps*span/2-1, sps), symbols, 'k.')
      plt.show()
    avg_pow = np.mean(np.abs(tx_stream)**2)
    tx_stream = tx_stream/np.sqrt(avg_pow)
    tx_stream_clean = tx_stream_clean/np.sqrt(avg_pow)
    symbols = symbols/np.sqrt(avg_pow)

    # need to make sure max sample value is less than 1 for USRP DAC input
    max_val = 1.01*np.max(tx_stream)
    tx_stream = tx_stream/max_val
    tx_stream_clean = tx_stream_clean/max_val
    symbols = symbols/max_val
    # import pdb; pdb.set_trace()

    return symbols, tx_stream, tx_stream_clean

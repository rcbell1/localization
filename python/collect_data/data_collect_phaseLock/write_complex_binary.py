import numpy as np

def write_complex_binary(samples, file_name):

    samples = np.array(samples)
    rsamples = samples.real.reshape((-1,1)).astype(np.float32)
    isamples = samples.imag.reshape((-1,1)).astype(np.float32)
    y = np.concatenate( (rsamples, isamples), axis=1 )
    y.tofile(file_name) # write binary data to file in row major order

if __name__ == '__main__':
    # these lines are for debug from command line
    samps = np.random.rand(10,1)+1j*np.random.rand(10,1)
    fname = 'test.dat'
    print(samps)
    write_complex_binary(samps, fname)

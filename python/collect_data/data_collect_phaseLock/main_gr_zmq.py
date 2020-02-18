# import zmq
import os
import sys
import numpy as np
import time
import usb.core
import usb.util
import usb.backend.libusb1 as libusb1
from tx_1n import tx_1n
from rx_3n_phaselock_pub import rx_3n_phaselock_pub
from rx_3n_phaselock_sub import rx_3n_phaselock_sub
from write_complex_binary import write_complex_binary
from tx_waveform import transmit_waveform

import signal
import pdb

#SIGINT handler
def sigint_handler(signal, frame):
    #Do something while breaking
    print("SIGINT called, exiting....")
    sys.exit(0)
#---------------------------------------------------------------------------

def setup_usb_switch(num_switches, default_port):
    # find our device
    # I found this backend
    # fix here: https://github.com/pyusb/pyusb/blob/master/docs/tutorial.rst
    backend = usb.backend.libusb1.get_backend(find_library=lambda x: "/usr/lib/x86_64-linux-gnu/libusb-1.0.so")
    dev = usb.core.find(find_all=True,idVendor=0x20ce, idProduct=0x0022, backend=backend)
    dev = list(dev) # dev is not subscriptable without this line

    #was it found?
    if dev is None:
        raise ValueError('Device not found')

    # Detach kernal driver
    for switch in range(num_switches):
        for configuration in dev[switch]:
            for interface in configuration:
                ifnum = interface.bInterfaceNumber
                if not dev[switch].is_kernel_driver_active(ifnum):
                    continue
                try:
                    dev[switch].detach_kernel_driver(ifnum)
                except usb.core.USBError:
                    pass

        dev[switch].set_configuration()
        dev[switch].write(1,chr(default_port)) # default set upon startup
        print("Switch {0} set to port {1} by default".format(switch, default_port))

    return dev

def main():
    signal.signal(signal.SIGINT, sigint_handler)
    npulses = 100
    num_zeros = 3000   # number of zeros between pulses in symbols
    num_prn = 1000  # number of prn symbols per pulse in symbols
    roll_off = 0.5  # the shaping filter roll off factor
    span = 10       # the span of shaping filter in symbols
    rx_samp_rate = 200e6/22
    tx_samp_rate = 200e6/66
    center_freq = 2.395e9
    sps = 2
    tx_gain = 25
    rx_gain = 18
    num_rx = 3
    socket_addr = 'tcp://localhost:8000'
    filepath_save = '~/Documents/repos/localization/data/14/tx_center/rfs9/'
    num_switches = 3
    default_port = 2
    zmq_pub_high_water_mark = 5
    zmq_sub_high_water_mark = 10
    flush_zmq_buffers_en = 1
    tx_enable = 1

    dev = setup_usb_switch(num_switches, default_port)
    tx_symbols, tx_pulse, tx_pulse_clean = \
        transmit_waveform(Nsym=num_prn,sps=sps,span=span,rolloff=roll_off,\
            pad_zeros=num_zeros,plot_debug=False)
    # tx_symbols = list(np.zeros(int(num_zeros/2)))+ \
    #     list(2*np.random.randint(0,2,size=num_prn)-1)+ \
    #     list(np.zeros(int(num_zeros/2)))

    # instantiate the GNU Radio flowgraph
    if tx_enable == 1:
        tx_stream = list(np.zeros(sps*int(num_zeros/2)))+ \
            list(tx_pulse) + \
            list(np.zeros(sps*int(num_zeros/2)))
        tx = tx_1n()
        tx.set_samp_rate(tx_samp_rate)
        tx.set_center_freq(center_freq)
        tx.set_sps(sps)
        tx.set_tx_gain(tx_gain)
        tx.set_num_zeros(num_zeros)
        tx.set_num_prn(num_prn)
        tx.set_roll_off(roll_off)
        tx.set_span(span)
        tx.blocks_vector_source_x_0.set_data(tx_stream)

    rx_pub = rx_3n_phaselock_pub()
    rx_pub.set_samp_rate(rx_samp_rate)
    rx_pub.set_center_freq(center_freq)
    rx_pub.set_rx_gain(rx_gain)
    rx_pub.set_socket_addr(socket_addr)
    rx_pub.set_high_water_mark(zmq_pub_high_water_mark)

    rx_sub = rx_3n_phaselock_sub()
    rx_sub.set_npulses_stop(npulses)
    rx_sub.set_samp_rate(rx_samp_rate)
    rx_sub.set_socket_addr(socket_addr)
    rx_sub.set_tx_samp_rate(tx_samp_rate)
    rx_sub.set_sps(sps)
    rx_sub.set_high_water_mark(zmq_sub_high_water_mark)
    # print(rx_sub.get_high_water_mark())

    if tx_enable == 1:
        tx.start()
    time.sleep(0.5)

    rx_pub.start()
    time.sleep(0.5)

    print('''\nTo collect a data set press the \'c\' key.
To set switches to port X enter \'s X\', where X = 1,2,3,4.
To end this session press \'q\'.''')

    first_run = 1 # first time through loop flag
    set_num = 1 # data collection set number
    while(True):
        user_input = raw_input("\nEnter your command: ")
        tokens = user_input.split()

        if len(tokens) == 1:
            command = tokens[0]
        elif len(tokens) == 2:
            command = tokens[0]
            port_number = tokens[1]
        else:
            print('\nUnrecognized entry. Try Again')
            pass

        if command == 'c':
            print("\nCollecting {0} pulses\n".format(npulses))

            # ignore transients from startup
            if first_run == 1:
                first_run = 0

                rx_sub.start()
                rx_sub.wait()

                junk0 = rx_sub.blocks_vector_sink_x_0.data()
                junk1 = rx_sub.blocks_vector_sink_x_1.data()
                junk2 = rx_sub.blocks_vector_sink_x_2.data()
                rx_sub.blocks_head_0.reset()
                rx_sub.blocks_head_1.reset()
                rx_sub.blocks_head_2.reset()
                rx_sub.blocks_vector_sink_x_0.reset()
                rx_sub.blocks_vector_sink_x_1.reset()
                rx_sub.blocks_vector_sink_x_2.reset()

            if flush_zmq_buffers_en == 1:
                for flush in range(10):
                    rx_sub.start()
                    rx_sub.wait()

                    junk0 = rx_sub.blocks_vector_sink_x_0.data()
                    junk1 = rx_sub.blocks_vector_sink_x_1.data()
                    junk2 = rx_sub.blocks_vector_sink_x_2.data()
                    rx_sub.blocks_head_0.reset()
                    rx_sub.blocks_head_1.reset()
                    rx_sub.blocks_head_2.reset()
                    rx_sub.blocks_vector_sink_x_0.reset()
                    rx_sub.blocks_vector_sink_x_1.reset()
                    rx_sub.blocks_vector_sink_x_2.reset()


            rx_sub.start()
            rx_sub.wait()

            # get the data from vector sinks in GNU Radio
            data0 = rx_sub.blocks_vector_sink_x_0.data()
            data1 = rx_sub.blocks_vector_sink_x_1.data()
            data2 = rx_sub.blocks_vector_sink_x_2.data()

            # save the collected data to files
            path = filepath_save + "{0}".format(set_num) + '/'
            try:
                os_path = os.path.expanduser(path)
                os.makedirs(os_path)
            except OSError:
                if not os.path.isdir(os_path):
                    raise

            if tx_enable == 1:
                write_complex_binary(tx_pulse, os_path + 'tx_pulse.dat')
                write_complex_binary(tx_pulse_clean, os_path + 'tx_pulse_clean.dat')
                write_complex_binary(tx_symbols, os_path + 'tx_symbols.dat')
            write_complex_binary(data0, os_path + 'rx1.dat')
            write_complex_binary(data1, os_path + 'rx2.dat')
            write_complex_binary(data2, os_path + 'rx3.dat')
            set_num += 1

            # reset the head blocks for next collection
            rx_sub.blocks_head_0.reset()
            rx_sub.blocks_head_1.reset()
            rx_sub.blocks_head_2.reset()
            rx_sub.blocks_vector_sink_x_0.reset()
            rx_sub.blocks_vector_sink_x_1.reset()
            rx_sub.blocks_vector_sink_x_2.reset()

            print("Rx0 Samples:\n {0}\n".format(data0[0:9]))
            print("Rx1 Samples:\n {0}\n".format(data1[0:9]))
            print("Rx2 Samples:\n {0}\n".format(data2[0:9]))
            print("Files saved to: {0}\n\n".format(os_path))

            # rx_sub.stop()
            # rx_sub.wait()

        elif command == 's':
            if port_number in ['1','2','3','4']:
                for switch in range(num_switches):
                    dev[switch].write(1,chr(int(port_number)))
                    print("Switch {0} to port {1}".format(switch, port_number))
            else:
                print('\nUnrecognized entry. Try Again')
                pass

        elif command == 'q':
            if tx_enable == 1:
                tx.stop()
                tx.wait()
            rx_pub.stop()
            rx_pub.wait()
            rx_sub.stop()
            rx_sub.wait()
            print('exit')
            # exit(0)
            break

        else:
            print('\nUnrecognized entry. Try Again')
            pass

    print "\nDone\n"

if __name__ == '__main__':
    main()

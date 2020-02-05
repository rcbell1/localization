# import zmq
import os
import numpy as np
import matplotlib.pyplot as plt
import time
import usb.core
import usb.util
import usb.backend.libusb1 as libusb1
from rx_3n_phaselock_pub import rx_3n_phaselock_pub
from rx_3n_phaselock_sub import rx_3n_phaselock_sub
from write_complex_binary import write_complex_binary
# from switch_sp4t_63 import setup_usb_switch

import pdb

def setup_usb_switch(num_switches):
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

    return dev

def main():
    npulses = 100
    samp_rate = 200e6/22
    center_freq = 2.395e9
    rx_gain = 25
    num_rx = 3
    socket_addr = 'tcp://localhost:8000'
    filepath_save = '~/Documents/repos/localization/data/13/tx_center/rfs9/'
    num_switches = 3

    dev = setup_usb_switch(num_switches)

    # instantiate the GNU Radio flowgraph
    rx_pub = rx_3n_phaselock_pub()
    rx_pub.set_samp_rate(samp_rate)
    rx_pub.set_center_freq(center_freq)
    rx_pub.set_rx_gain(rx_gain)
    rx_pub.set_socket_addr(socket_addr)

    rx_sub = rx_3n_phaselock_sub()
    rx_sub.set_npulses_stop(npulses)
    rx_sub.set_samp_rate(samp_rate)
    rx_sub.set_socket_addr(socket_addr)

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
            switch_number = tokens[1]
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
            print("Rx2 Samples:\n {0}\n\n".format(data2[0:9]))

        elif command == 's':
            if switch_number in ['1','2','3','4']:
                dev[0].write(1,switch_number) #v11=1 switch port 1; v11=2 switch port 2
                dev[1].write(1,switch_number) #v22=1 switch port 1; v22=2 switch port 2
                dev[2].write(1,switch_number)
                print("Switched to port {0}".format(switch_number))
            else:
                print('\nUnrecognized entry. Try Again')
                pass

        elif command == 'q':
            rx_pub.stop()
            rx_pub.wait()
            rx_sub.stop()
            rx_sub.wait()
            break

        else:
            print('\nUnrecognized entry. Try Again')
            pass

    print "\nDone\n"

if __name__ == '__main__':
    main()

import zmq
import numpy as np
import matplotlib.pyplot as plt
import time
import usb.core
import usb.util
import usb.backend.libusb1 as libusb1
from rx_3n import rx_3n
from write_complex_binary import write_complex_binary

import pdb

def setup_usb_switch():
    # find our device
    # I found this backend
    # fix here: https://github.com/pyusb/pyusb/blob/master/docs/tutorial.rst
    backend = usb.backend.libusb1.get_backend(find_library=lambda x: "/usr/lib/x86_64-linux-gnu/libusb-1.0.so")
    dev = usb.core.find(find_all=True,idVendor=0x20ce, idProduct=0x0022, backend=backend)
    dev = list(dev) # dev is not subscriptable without this line

    #dev = usb.core.find(idVendor=0x20ce, idProduct=0x0022)
    # langid = 0x0000
    # breakpoint()
    # print(usb.util.get_string(dev[0], 256, dev[0].iSerialNumber))
    # print(usb.util.get_string(dev[1], 256, dev[1].iSerialNumber))

    #was it found?
    if dev is None:
        raise ValueError('Device not found')

    # Detach kernal driver
    for configuration in dev[0]:
        for interface in configuration:
            ifnum = interface.bInterfaceNumber
            if not dev[0].is_kernel_driver_active(ifnum):
                continue
            try:
                dev[0].detach_kernel_driver(ifnum)
            except usb.core.USBError:
                pass

    for configuration in dev[1]:
        for interface in configuration:
            ifnum = interface.bInterfaceNumber
            if not dev[1].is_kernel_driver_active(ifnum):
                continue
            try:
                dev[1].detach_kernel_driver(ifnum)
            except usb.core.USBError:
                pass

    #set the active configuration. with no args we use first configure
    dev[0].set_configuration()
    dev[1].set_configuration()

    return dev

def main():
    npulses = 100
    samp_rate = 200e6/22
    center_freq = 2.395e9
    rx_gain = 25
    num_rx = 3

    dev = setup_usb_switch()

    # instantiate the GNU Radio flowgraph
    receiver = rx_3n()
    receiver.set_npulses_stop(npulses)
    receiver.set_samp_rate(samp_rate)
    receiver.set_center_freq(center_freq)
    receiver.set_rx_gain(rx_gain)
    time.sleep(1)

    print('''\nTo collect a data set press the \'c\' key.
To set switches to port X enter \'s X\', where X = 1,2,3,4.
To end this session press \'q\'.''')

    first_run = 1 # first time through loop flag
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
                receiver.start()
                receiver.wait()
                junk0 = receiver.blocks_vector_sink_x_0.data()
                junk1 = receiver.blocks_vector_sink_x_1.data()
                junk2 = receiver.blocks_vector_sink_x_2.data()
                receiver.blocks_head_0.reset()
                receiver.blocks_head_1.reset()
                receiver.blocks_head_2.reset()
                receiver.blocks_vector_sink_x_0.reset()
                receiver.blocks_vector_sink_x_1.reset()
                receiver.blocks_vector_sink_x_2.reset()

            receiver.start()
            receiver.wait()

            # get the data from vector sinks in GNU Radio
            data0 = receiver.blocks_vector_sink_x_0.data()
            data1 = receiver.blocks_vector_sink_x_1.data()
            data2 = receiver.blocks_vector_sink_x_2.data()

            # save the collected data to files
            write_complex_binary(data0, 'rx1.dat')
            write_complex_binary(data0, 'rx2.dat')
            write_complex_binary(data0, 'rx3.dat')

            # reset the head blocks for next collection
            receiver.blocks_head_0.reset()
            receiver.blocks_head_1.reset()
            receiver.blocks_head_2.reset()
            receiver.blocks_vector_sink_x_0.reset()
            receiver.blocks_vector_sink_x_1.reset()
            receiver.blocks_vector_sink_x_2.reset()

            print("Rx0 Samples:\n {0}\n".format(data0[0:9]))
            print("Rx1 Samples:\n {0}\n".format(data1[0:9]))
            print("Rx2 Samples:\n {0}\n\n".format(data2[0:9]))

        elif command == 's':
            if switch_number in ['1','2','3','4']:
                dev[0].write(1,switch_number) #v11=1 switch port 1; v11=2 switch port 2
                dev[1].write(1,switch_number) #v22=1 switch port 1; v22=2 switch port 2
                print("Switched to port {0}".format(switch_number))
            else:
                print('\nUnrecognized entry. Try Again')
                pass

        elif command == 'q':
            break

        else:
            print('\nUnrecognized entry. Try Again')
            pass

    receiver.stop()
    print "\nDone\n"

if __name__ == '__main__':
    main()

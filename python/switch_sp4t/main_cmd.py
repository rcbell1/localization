import usb.core
import usb.util
import usb.backend.libusb1 as libusb1

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
        print("Switch {0} set to port {1} by default".format(switch,'2'))

    return dev

def main():
    num_switches = 3
    default_port = 1

    dev = setup_usb_switch(num_switches, default_port)

    print('''\nTo set switches to port X enter \'s X\', where X = 1,2,3,4.
To end this session press \'q\'.''')

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

        if command == 's':
            if port_number in ['1','2','3','4']:
                for switch in range(num_switches):
                    dev[switch].write(1,chr(int(port_number))) 
                    print("Switch {0} to port {1}".format(switch, port_number))
            else:
                print('\nUnrecognized entry. Try Again')
                pass

        elif command == 'q':
            break

        else:
            print('\nUnrecognized entry. Try Again')
            pass

    print "\nDone\n"

if __name__ == '__main__':
    main()

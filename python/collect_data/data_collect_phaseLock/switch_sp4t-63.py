import usb.core
import usb.util
import usb.backend.libusb1 as libusb1

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

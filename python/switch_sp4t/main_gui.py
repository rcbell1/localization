# This code was pulled from section 2.1.3 (e) - Solid State Switches - 2 Devices
# found here:
# https://www.minicircuits.com/softwaredownload/Prog_Examples_Troubleshooting.pdf

import time
import usb.core
import usb.util
import tkinter as tk
import usb.backend.libusb1 as libusb1

def ApplySwitching():
    v11 = int(V1.get())
    v22 = int(V2.get())
    dev[0].write(1,chr(v11)) #v11=1 switch port 1; v11=2 switch port 2
    dev[1].write(1,chr(v22)) #v22=1 switch port 1; v22=2 switch port 2
    # import pdb; pdb.set_trace()
    print(v11, v22)

root = tk.Tk()
root.title("Mini Circuits USB-SP4T ver 1.0")
tk.Label(root, text="    ").grid(row=2, column=1)
V1 = tk.IntVar()
V2 = tk.IntVar()

# breakpoint()

tk.Button(root, text = "Apply Switching", width=10,height=1,command=ApplySwitching).grid(row=4,column=1)
tk.Spinbox(root,values=(1,2,3,4),width=5,textvariable=V1).grid(row=3,column=2)
tk.Spinbox(root,values=(1,2,3,4),width=5,textvariable=V2).grid(row=3,column=3)

tk.Label(root,text="Comm -> Switch: ").grid(row=3,column=1)

#find our device
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

#set the active configuration. with no args we usue first configure
dev[0].set_configuration()
dev[1].set_configuration()

root.geometry("500x300+0+0")
root.mainloop()

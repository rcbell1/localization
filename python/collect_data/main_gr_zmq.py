import zmq
import numpy as np
import matplotlib.pyplot as plt
import time
from rx_3n import rx_3n
from tb_zero_mq import tb_zero_mq

import pdb


def main():
    ip_addr_0 = "localhost"; port0 = "8000"
    ip_addr_1 = "localhost"; port1 = "8001"
    ip_addr_2 = "localhost"; port2 = "8002"

    # instantiate the GNU Radio flowgraph which will act as a zmq data server
    # receiver = rx_3n()
    # receiver.set_socket_addr_0("tcp://{0}:{1}".format(ip_addr_0, port0))
    # receiver.set_socket_addr_1("tcp://{0}:{1}".format(ip_addr_1, port1))
    # receiver.set_socket_addr_2("tcp://{0}:{1}".format(ip_addr_2, port2))

    tb = tb_zero_mq()
    tb.set_socket_addr("tcp://{0}:{1}".format(ip_addr_0, port0))

    # create the zmq client that requests data from the GNU Radio server
    context = zmq.Context()
    print "Connecting to server..."
    socket = context.socket(zmq.REQ)
    socket.connect ("tcp://{0}:{1}".format(ip_addr_0, port0))
    # socket.connect ("tcp://{0}:{1}".format(ip_addr_1, port1))
    # socket.connect ("tcp://{0}:{1}".format(ip_addr_2, port2))

    # receiver.start()
    tb.start()
    time.sleep(2)

    for request in range(1,2):
        print "Sending request ", request,"..."
        socket.send ("Hello")
        #  Get the reply.
        message = socket.recv()
        pdb.set_trace()
        print "Received reply ", request, "[", message, "]"

    receiver.stop()
    # zmq.close(socket)
    socket.close()
    context.destroy()
    print "Done"

if __name__ == '__main__':
    main()

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Rx 3N Phaselock Sub
# GNU Radio version: 3.7.13.5
##################################################

from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import zeromq
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import numpy as np


class rx_3n_phaselock_sub(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Rx 3N Phaselock Sub")

        ##################################################
        # Variables
        ##################################################
        self.tx_samp_rate = tx_samp_rate = 200e6/70
        self.sps = sps = 4
        self.samp_rate = samp_rate = 200e6/12
        self.prnLen = prnLen = 1000
        self.nzeros = nzeros = 3000
        self.npulses_stop = npulses_stop = 100
        self.span = span = 10
        self.socket_addr = socket_addr = 'tcp://127.0.0.1:8000'
        self.nitems_stop = nitems_stop = np.ceil(sps*(prnLen+nzeros)*npulses_stop*samp_rate/tx_samp_rate)

        ##################################################
        # Blocks
        ##################################################
        self.zeromq_sub_source_0 = zeromq.sub_source(gr.sizeof_gr_complex, 1, 'tcp://127.0.0.1:6000', 100, False, -1)
        self.blocks_vector_sink_x_2 = blocks.vector_sink_c(1, 1024)
        self.blocks_vector_sink_x_1 = blocks.vector_sink_c(1, 1024)
        self.blocks_vector_sink_x_0 = blocks.vector_sink_c(1, 1024)
        self.blocks_stream_to_streams_0 = blocks.stream_to_streams(gr.sizeof_gr_complex*1, 3)
        self.blocks_head_2 = blocks.head(gr.sizeof_gr_complex*1, int(nitems_stop))
        self.blocks_head_1 = blocks.head(gr.sizeof_gr_complex*1, int(nitems_stop))
        self.blocks_head_0 = blocks.head(gr.sizeof_gr_complex*1, int(nitems_stop))



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_head_0, 0), (self.blocks_vector_sink_x_0, 0))
        self.connect((self.blocks_head_1, 0), (self.blocks_vector_sink_x_1, 0))
        self.connect((self.blocks_head_2, 0), (self.blocks_vector_sink_x_2, 0))
        self.connect((self.blocks_stream_to_streams_0, 0), (self.blocks_head_0, 0))
        self.connect((self.blocks_stream_to_streams_0, 1), (self.blocks_head_1, 0))
        self.connect((self.blocks_stream_to_streams_0, 2), (self.blocks_head_2, 0))
        self.connect((self.zeromq_sub_source_0, 0), (self.blocks_stream_to_streams_0, 0))

    def get_tx_samp_rate(self):
        return self.tx_samp_rate

    def set_tx_samp_rate(self, tx_samp_rate):
        self.tx_samp_rate = tx_samp_rate
        self.set_nitems_stop(np.ceil(self.sps*(self.prnLen+self.nzeros)*self.npulses_stop*self.samp_rate/self.tx_samp_rate))

    def get_sps(self):
        return self.sps

    def set_sps(self, sps):
        self.sps = sps
        self.set_nitems_stop(np.ceil(self.sps*(self.prnLen+self.nzeros)*self.npulses_stop*self.samp_rate/self.tx_samp_rate))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.set_nitems_stop(np.ceil(self.sps*(self.prnLen+self.nzeros)*self.npulses_stop*self.samp_rate/self.tx_samp_rate))

    def get_prnLen(self):
        return self.prnLen

    def set_prnLen(self, prnLen):
        self.prnLen = prnLen
        self.set_nitems_stop(np.ceil(self.sps*(self.prnLen+self.nzeros)*self.npulses_stop*self.samp_rate/self.tx_samp_rate))

    def get_nzeros(self):
        return self.nzeros

    def set_nzeros(self, nzeros):
        self.nzeros = nzeros
        self.set_nitems_stop(np.ceil(self.sps*(self.prnLen+self.nzeros)*self.npulses_stop*self.samp_rate/self.tx_samp_rate))

    def get_npulses_stop(self):
        return self.npulses_stop

    def set_npulses_stop(self, npulses_stop):
        self.npulses_stop = npulses_stop
        self.set_nitems_stop(np.ceil(self.sps*(self.prnLen+self.nzeros)*self.npulses_stop*self.samp_rate/self.tx_samp_rate))

    def get_span(self):
        return self.span

    def set_span(self, span):
        self.span = span

    def get_socket_addr(self):
        return self.socket_addr

    def set_socket_addr(self, socket_addr):
        self.socket_addr = socket_addr

    def get_nitems_stop(self):
        return self.nitems_stop

    def set_nitems_stop(self, nitems_stop):
        self.nitems_stop = nitems_stop
        self.blocks_head_2.set_length(int(self.nitems_stop))
        self.blocks_head_1.set_length(int(self.nitems_stop))
        self.blocks_head_0.set_length(int(self.nitems_stop))


def main(top_block_cls=rx_3n_phaselock_sub, options=None):

    tb = top_block_cls()
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Tx 1N
# GNU Radio version: 3.7.13.5
##################################################

from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import numpy as np
import time


class tx_1n(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Tx 1N")

        ##################################################
        # Variables
        ##################################################
        self.tx_gain = tx_gain = 0
        self.sps = sps = 4
        self.span = span = 20
        self.samp_rate = samp_rate = 2e6
        self.prnLen = prnLen = 1000
        self.nzeros = nzeros = 1
        self.center_freq = center_freq = 500e6

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_sink_0 = uhd.usrp_sink(
        	",".join(("addr0=192.168.10.3", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_sink_0.set_samp_rate(samp_rate)
        self.uhd_usrp_sink_0.set_center_freq(center_freq, 0)
        self.uhd_usrp_sink_0.set_gain(tx_gain, 0)
        self.uhd_usrp_sink_0.set_antenna('TX/RX', 0)
        self.root_raised_cosine_filter_0 = filter.interp_fir_filter_ccf(sps, firdes.root_raised_cosine(
        	0.8, sps, 1, 0.35, sps*span))
        self.blocks_vector_source_x_0 = blocks.vector_source_c(list(np.zeros(nzeros))+list(2*np.random.randint(0,2,size=prnLen)-1)+list(np.zeros(nzeros)), True, 1, [])



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_vector_source_x_0, 0), (self.root_raised_cosine_filter_0, 0))
        self.connect((self.root_raised_cosine_filter_0, 0), (self.uhd_usrp_sink_0, 0))

    def get_tx_gain(self):
        return self.tx_gain

    def set_tx_gain(self, tx_gain):
        self.tx_gain = tx_gain
        self.uhd_usrp_sink_0.set_gain(self.tx_gain, 0)


    def get_sps(self):
        return self.sps

    def set_sps(self, sps):
        self.sps = sps
        self.root_raised_cosine_filter_0.set_taps(firdes.root_raised_cosine(0.8, self.sps, 1, 0.35, self.sps*self.span))

    def get_span(self):
        return self.span

    def set_span(self, span):
        self.span = span
        self.root_raised_cosine_filter_0.set_taps(firdes.root_raised_cosine(0.8, self.sps, 1, 0.35, self.sps*self.span))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_sink_0.set_samp_rate(self.samp_rate)

    def get_prnLen(self):
        return self.prnLen

    def set_prnLen(self, prnLen):
        self.prnLen = prnLen
        self.blocks_vector_source_x_0.set_data(list(np.zeros(self.nzeros))+list(2*np.random.randint(0,2,size=self.prnLen)-1)+list(np.zeros(self.nzeros)), [])

    def get_nzeros(self):
        return self.nzeros

    def set_nzeros(self, nzeros):
        self.nzeros = nzeros
        self.blocks_vector_source_x_0.set_data(list(np.zeros(self.nzeros))+list(2*np.random.randint(0,2,size=self.prnLen)-1)+list(np.zeros(self.nzeros)), [])

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.uhd_usrp_sink_0.set_center_freq(self.center_freq, 0)


def main(top_block_cls=tx_1n, options=None):

    tb = top_block_cls()
    tb.start()
    try:
        raw_input('Press Enter to quit: ')
    except EOFError:
        pass
    tb.stop()
    tb.wait()


if __name__ == '__main__':
    main()

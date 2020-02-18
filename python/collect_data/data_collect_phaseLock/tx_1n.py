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
        self.tx_gain = tx_gain = 25
        self.sps = sps = 2
        self.span = span = 10
        self.samp_rate = samp_rate = 200e6/70
        self.roll_off = roll_off = 0.5
        self.num_zeros = num_zeros = 3000
        self.num_prn = num_prn = 1000
        self.center_freq = center_freq = 2.395e9

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_sink_0 = uhd.usrp_sink(
        	",".join(("addr0=192.168.13.2", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_sink_0.set_samp_rate(samp_rate)
        self.uhd_usrp_sink_0.set_center_freq(center_freq, 0)
        self.uhd_usrp_sink_0.set_gain(tx_gain, 0)
        self.uhd_usrp_sink_0.set_antenna('TX/RX', 0)
        self.blocks_vector_source_x_0 = blocks.vector_source_c(list(np.zeros(int(num_zeros/2)))+list(2*np.random.randint(0,2,size=num_prn)-1)+list(np.zeros(int(num_zeros/2))), True, 1, [])



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_vector_source_x_0, 0), (self.uhd_usrp_sink_0, 0))

    def get_tx_gain(self):
        return self.tx_gain

    def set_tx_gain(self, tx_gain):
        self.tx_gain = tx_gain
        self.uhd_usrp_sink_0.set_gain(self.tx_gain, 0)


    def get_sps(self):
        return self.sps

    def set_sps(self, sps):
        self.sps = sps

    def get_span(self):
        return self.span

    def set_span(self, span):
        self.span = span

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_sink_0.set_samp_rate(self.samp_rate)

    def get_roll_off(self):
        return self.roll_off

    def set_roll_off(self, roll_off):
        self.roll_off = roll_off

    def get_num_zeros(self):
        return self.num_zeros

    def set_num_zeros(self, num_zeros):
        self.num_zeros = num_zeros
        self.blocks_vector_source_x_0.set_data(list(np.zeros(int(self.num_zeros/2)))+list(2*np.random.randint(0,2,size=self.num_prn)-1)+list(np.zeros(int(self.num_zeros/2))), [])

    def get_num_prn(self):
        return self.num_prn

    def set_num_prn(self, num_prn):
        self.num_prn = num_prn
        self.blocks_vector_source_x_0.set_data(list(np.zeros(int(self.num_zeros/2)))+list(2*np.random.randint(0,2,size=self.num_prn)-1)+list(np.zeros(int(self.num_zeros/2))), [])

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.uhd_usrp_sink_0.set_center_freq(self.center_freq, 0)


def main(top_block_cls=tx_1n, options=None):

    tb = top_block_cls()
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Rx 3N
# Generated: Wed Dec 18 17:59:59 2019
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


class rx_3n(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Rx 3N")

        ##################################################
        # Variables
        ##################################################
        self.tx_samp_rate = tx_samp_rate = 200e6/70
        self.sps = sps = 4
        self.samp_rate = samp_rate = 200e6/40
        self.prnLen = prnLen = 1000
        self.nzeros = nzeros = 3000
        self.npulses_stop = npulses_stop = 100
        self.span = span = 10
        self.rx_gain = rx_gain = 30
        self.nitems_stop = nitems_stop = np.ceil(sps*(prnLen+nzeros)*npulses_stop*samp_rate/tx_samp_rate)
        self.center_freq = center_freq = 2.395e9

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_source_0 = uhd.usrp_source(
        	",".join(("addr0=192.168.10.2,addr1=192.168.10.4,addr2=192.168.10.7", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(3),
        	),
        )
        self.uhd_usrp_source_0.set_clock_source('external', 0)
        self.uhd_usrp_source_0.set_time_source('external', 0)
        self.uhd_usrp_source_0.set_subdev_spec('B:0', 0)
        self.uhd_usrp_source_0.set_clock_source('external', 1)
        self.uhd_usrp_source_0.set_time_source('external', 1)
        self.uhd_usrp_source_0.set_subdev_spec('B:0', 1)
        self.uhd_usrp_source_0.set_clock_source('external', 2)
        self.uhd_usrp_source_0.set_time_source('external', 2)
        self.uhd_usrp_source_0.set_subdev_spec('B:0', 2)
        self.uhd_usrp_source_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_source_0.set_center_freq(center_freq, 0)
        self.uhd_usrp_source_0.set_gain(rx_gain, 0)
        self.uhd_usrp_source_0.set_antenna('RX2', 0)
        self.uhd_usrp_source_0.set_center_freq(center_freq, 1)
        self.uhd_usrp_source_0.set_gain(rx_gain, 1)
        self.uhd_usrp_source_0.set_antenna('TX/RX', 1)
        self.uhd_usrp_source_0.set_center_freq(center_freq, 2)
        self.uhd_usrp_source_0.set_gain(rx_gain, 2)
        self.uhd_usrp_source_0.set_antenna('RX2', 2)
        self.dc_blocker_xx_0_0_0 = filter.dc_blocker_cc(200, True)
        self.dc_blocker_xx_0_0 = filter.dc_blocker_cc(200, True)
        self.dc_blocker_xx_0 = filter.dc_blocker_cc(200, True)
        self.blocks_head_0_1 = blocks.head(gr.sizeof_gr_complex*1, int(nitems_stop))
        self.blocks_head_0_0 = blocks.head(gr.sizeof_gr_complex*1, int(nitems_stop))
        self.blocks_head_0 = blocks.head(gr.sizeof_gr_complex*1, int(nitems_stop))
        self.blocks_file_sink_0_0_0 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/wcsng-21/Documents/richbell/tdoa-localization/data/2/tx_center/rx_samp_rate_5Msps/rx3.dat', False)
        self.blocks_file_sink_0_0_0.set_unbuffered(False)
        self.blocks_file_sink_0_0 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/wcsng-21/Documents/richbell/tdoa-localization/data/2/tx_center/rx_samp_rate_5Msps/rx2.dat', False)
        self.blocks_file_sink_0_0.set_unbuffered(False)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/wcsng-21/Documents/richbell/tdoa-localization/data/2/tx_center/rx_samp_rate_5Msps/rx1.dat', False)
        self.blocks_file_sink_0.set_unbuffered(False)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_head_0, 0), (self.blocks_file_sink_0, 0))
        self.connect((self.blocks_head_0_0, 0), (self.blocks_file_sink_0_0, 0))
        self.connect((self.blocks_head_0_1, 0), (self.blocks_file_sink_0_0_0, 0))
        self.connect((self.dc_blocker_xx_0, 0), (self.blocks_head_0, 0))
        self.connect((self.dc_blocker_xx_0_0, 0), (self.blocks_head_0_0, 0))
        self.connect((self.dc_blocker_xx_0_0_0, 0), (self.blocks_head_0_1, 0))
        self.connect((self.uhd_usrp_source_0, 0), (self.dc_blocker_xx_0, 0))
        self.connect((self.uhd_usrp_source_0, 1), (self.dc_blocker_xx_0_0, 0))
        self.connect((self.uhd_usrp_source_0, 2), (self.dc_blocker_xx_0_0_0, 0))

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
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)

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

    def get_rx_gain(self):
        return self.rx_gain

    def set_rx_gain(self, rx_gain):
        self.rx_gain = rx_gain
        self.uhd_usrp_source_0.set_gain(self.rx_gain, 0)

        self.uhd_usrp_source_0.set_gain(self.rx_gain, 1)

        self.uhd_usrp_source_0.set_gain(self.rx_gain, 2)


    def get_nitems_stop(self):
        return self.nitems_stop

    def set_nitems_stop(self, nitems_stop):
        self.nitems_stop = nitems_stop
        self.blocks_head_0_1.set_length(int(self.nitems_stop))
        self.blocks_head_0_0.set_length(int(self.nitems_stop))
        self.blocks_head_0.set_length(int(self.nitems_stop))

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.uhd_usrp_source_0.set_center_freq(self.center_freq, 0)
        self.uhd_usrp_source_0.set_center_freq(self.center_freq, 1)
        self.uhd_usrp_source_0.set_center_freq(self.center_freq, 2)


def main(top_block_cls=rx_3n, options=None):

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

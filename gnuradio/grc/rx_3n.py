#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Rx 3N
# Generated: Fri Nov 22 19:01:27 2019
##################################################

from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import time


class rx_3n(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Rx 3N")

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 200e6/22
        self.rx_gain = rx_gain = 30
        self.nitems_stop = nitems_stop = int(1e6)
        self.center_freq = center_freq = 2.5e9

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_source_0 = uhd.usrp_source(
        	",".join(("addr0=192.168.10.2,addr1=192.168.10.4,addr2=192.168.10.5", "")),
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
        self.blocks_head_0 = blocks.head(gr.sizeof_gr_complex*1, nitems_stop)
        self.blocks_file_sink_0_0_0 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/wcsng-21/Documents/richbell/tdoa-localization/data/tx_near_rx1/rx3.dat', False)
        self.blocks_file_sink_0_0_0.set_unbuffered(False)
        self.blocks_file_sink_0_0 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/wcsng-21/Documents/richbell/tdoa-localization/data/tx_near_rx1/rx1.dat', False)
        self.blocks_file_sink_0_0.set_unbuffered(False)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_gr_complex*1, '/home/wcsng-21/Documents/richbell/tdoa-localization/data/tx_near_rx1/rx2.dat', False)
        self.blocks_file_sink_0.set_unbuffered(False)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_head_0, 0), (self.blocks_file_sink_0, 0))
        self.connect((self.dc_blocker_xx_0, 0), (self.blocks_head_0, 0))
        self.connect((self.dc_blocker_xx_0_0, 0), (self.blocks_file_sink_0_0, 0))
        self.connect((self.dc_blocker_xx_0_0_0, 0), (self.blocks_file_sink_0_0_0, 0))
        self.connect((self.uhd_usrp_source_0, 0), (self.dc_blocker_xx_0, 0))
        self.connect((self.uhd_usrp_source_0, 1), (self.dc_blocker_xx_0_0, 0))
        self.connect((self.uhd_usrp_source_0, 2), (self.dc_blocker_xx_0_0_0, 0))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)

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
        self.blocks_head_0.set_length(self.nitems_stop)

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

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Rx 3N Phaselock Pub
# GNU Radio version: 3.7.13.5
##################################################

from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio import uhd
from gnuradio import zeromq
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import time


class rx_3n_phaselock_pub(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Rx 3N Phaselock Pub")

        ##################################################
        # Variables
        ##################################################
        self.socket_addr = socket_addr = 'tcp://127.0.0.1:8000'
        self.samp_rate = samp_rate = 200e6/12
        self.rx_gain = rx_gain = 20
        self.high_water_mark = high_water_mark = -1
        self.center_freq = center_freq = 2.395e9

        ##################################################
        # Blocks
        ##################################################
        self.zeromq_pub_sink_0 = zeromq.pub_sink(gr.sizeof_gr_complex, 1, socket_addr, 100, False, high_water_mark)
        self.uhd_usrp_source_0 = uhd.usrp_source(
        	",".join(("addr0=192.168.10.2,addr1=192.168.11.2,addr2=192.168.12.2", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(3),
        	),
        )
        self.uhd_usrp_source_0.set_clock_rate(200e6, uhd.ALL_MBOARDS)
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
        self.uhd_usrp_source_0.set_auto_dc_offset(False, 0)
        self.uhd_usrp_source_0.set_auto_iq_balance(False, 0)
        self.uhd_usrp_source_0.set_center_freq(center_freq, 1)
        self.uhd_usrp_source_0.set_gain(rx_gain, 1)
        self.uhd_usrp_source_0.set_antenna('TX/RX', 1)
        self.uhd_usrp_source_0.set_auto_dc_offset(False, 1)
        self.uhd_usrp_source_0.set_auto_iq_balance(False, 1)
        self.uhd_usrp_source_0.set_center_freq(center_freq, 2)
        self.uhd_usrp_source_0.set_gain(rx_gain, 2)
        self.uhd_usrp_source_0.set_antenna('RX2', 2)
        self.uhd_usrp_source_0.set_auto_dc_offset(False, 2)
        self.uhd_usrp_source_0.set_auto_iq_balance(False, 2)
        self.dc_blocker_xx_0_0_0 = filter.dc_blocker_cc(200, True)
        self.dc_blocker_xx_0_0 = filter.dc_blocker_cc(200, True)
        self.dc_blocker_xx_0 = filter.dc_blocker_cc(200, True)
        self.blocks_streams_to_stream_0 = blocks.streams_to_stream(gr.sizeof_gr_complex*1, 3)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_streams_to_stream_0, 0), (self.zeromq_pub_sink_0, 0))
        self.connect((self.dc_blocker_xx_0, 0), (self.blocks_streams_to_stream_0, 0))
        self.connect((self.dc_blocker_xx_0_0, 0), (self.blocks_streams_to_stream_0, 1))
        self.connect((self.dc_blocker_xx_0_0_0, 0), (self.blocks_streams_to_stream_0, 2))
        self.connect((self.uhd_usrp_source_0, 0), (self.dc_blocker_xx_0, 0))
        self.connect((self.uhd_usrp_source_0, 1), (self.dc_blocker_xx_0_0, 0))
        self.connect((self.uhd_usrp_source_0, 2), (self.dc_blocker_xx_0_0_0, 0))

    def get_socket_addr(self):
        return self.socket_addr

    def set_socket_addr(self, socket_addr):
        self.socket_addr = socket_addr

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


    def get_high_water_mark(self):
        return self.high_water_mark

    def set_high_water_mark(self, high_water_mark):
        self.high_water_mark = high_water_mark

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.uhd_usrp_source_0.set_center_freq(self.center_freq, 0)
        self.uhd_usrp_source_0.set_center_freq(self.center_freq, 1)
        self.uhd_usrp_source_0.set_center_freq(self.center_freq, 2)


def main(top_block_cls=rx_3n_phaselock_pub, options=None):

    tb = top_block_cls()
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()

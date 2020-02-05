#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Tb Zero Mq
# GNU Radio version: 3.7.13.5
##################################################

from gnuradio import analog
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import zeromq
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser


class tb_zero_mq(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Tb Zero Mq")

        ##################################################
        # Variables
        ##################################################
        self.socket_addr = socket_addr = 'tcp://127.0.0.1:8000'
        self.samp_rate = samp_rate = 100

        ##################################################
        # Blocks
        ##################################################
        self.zeromq_rep_sink_0 = zeromq.rep_sink(gr.sizeof_short, 1, socket_addr, 100, False, 10)
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_short*1, samp_rate,True)
        self.analog_sig_source_x_0 = analog.sig_source_s(100, analog.GR_SAW_WAVE, 1, 100, 1)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_sig_source_x_0, 0), (self.blocks_throttle_0, 0))
        self.connect((self.blocks_throttle_0, 0), (self.zeromq_rep_sink_0, 0))

    def get_socket_addr(self):
        return self.socket_addr

    def set_socket_addr(self, socket_addr):
        self.socket_addr = socket_addr

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.blocks_throttle_0.set_sample_rate(self.samp_rate)


def main(top_block_cls=tb_zero_mq, options=None):

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

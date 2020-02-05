#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import time
import pdb

import csv

def main():

    samp_rate = 1e6;
    f_c = 4.35e9;
    rx_gain = 25;
    ant = 'RX2'
    # I use a B210 so you should include the correct IP
    # use: join(("addr",""))
    addr = "addr=192.168.10.1"

    uhd_usrp_source_0 = uhd.usrp_source(
            ",".join(("", "")),
            uhd.stream_args(
                cpu_format="fc32",
                channels=range(1),
            ),
        )
    uhd_usrp_source_0.set_samp_rate(samp_rate)
    uhd_usrp_source_0.set_center_freq(f_c, 0)
    uhd_usrp_source_0.set_gain(rx_gain, 0)
    uhd_usrp_source_0.set_antenna(ant, 0)
    uhd_usrp_source_0.set_auto_dc_offset(True, 0)
    uhd_usrp_source_0.set_auto_iq_balance(True, 0)

    ax1 = uhd_usrp_source_0.finite_acquisition(1000000)

    ##
    # Put switching code here
    #
    ##

    ax2 = uhd_usrp_source_0.finite_acquisition(1000000)

    with open('filename.csv','wb') as out:
        csv_out=csv.writer(out)
        for samp in ax1:
            csv_out.writerow([samp.real,samp.imag])
        for samp in ax2:
            csv_out.writerow([samp.real,samp.imag])
    # pdb.set_trace()


if __name__ == '__main__':
    main()

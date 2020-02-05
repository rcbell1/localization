#!/usr/bin/env python3
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
    num_samps_collect = 10000
    samp_rate = 200e6/22
    center_freq = 2.395e9
    rx_gain = 25
    # I use a B210 so you should include the correct IP
    # use: join(("addr",""))
#     addr = "addr=192.168.10.1"
# ​
#     uhd_usrp_source_0 = uhd.usrp_source(
#             ",".join(("", "")),
#             uhd.stream_args(
#                 cpu_format="fc32",
#                 channels=range(1),
#             ),
#         )
#     uhd_usrp_source_0.set_samp_rate(samp_rate)
#     uhd_usrp_source_0.set_center_freq(f_c, 0)
#     uhd_usrp_source_0.set_gain(rx_gain, 0)
#     uhd_usrp_source_0.set_antenna(ant, 0)
#     uhd_usrp_source_0.set_auto_dc_offset(True, 0)
#     uhd_usrp_source_0.set_auto_iq_balance(True, 0)
    uhd_usrp_source_0 = uhd.usrp_source(
    	",".join(("addr0=192.168.10.2,addr1=192.168.11.2,addr2=192.168.12.2", "")),
    	uhd.stream_args(
    		cpu_format="fc32",
    		channels=range(3),
    	),
    )
    uhd_usrp_source_0.set_clock_rate(200e6, uhd.ALL_MBOARDS)
    uhd_usrp_source_0.set_clock_source('external', 0)
    uhd_usrp_source_0.set_time_source('external', 0)
    uhd_usrp_source_0.set_subdev_spec('B:0', 0)
    uhd_usrp_source_0.set_clock_source('external', 1)
    uhd_usrp_source_0.set_time_source('external', 1)
    uhd_usrp_source_0.set_subdev_spec('B:0', 1)
    uhd_usrp_source_0.set_clock_source('external', 2)
    uhd_usrp_source_0.set_time_source('external', 2)
    uhd_usrp_source_0.set_subdev_spec('B:0', 2)
    uhd_usrp_source_0.set_samp_rate(samp_rate)
    uhd_usrp_source_0.set_time_unknown_pps(uhd.time_spec())
    uhd_usrp_source_0.set_center_freq(center_freq, 0)
    uhd_usrp_source_0.set_gain(rx_gain, 0)
    uhd_usrp_source_0.set_antenna('RX2', 0)
    uhd_usrp_source_0.set_auto_dc_offset(False, 0)
    uhd_usrp_source_0.set_auto_iq_balance(False, 0)
    uhd_usrp_source_0.set_center_freq(center_freq, 1)
    uhd_usrp_source_0.set_gain(rx_gain, 1)
    uhd_usrp_source_0.set_antenna('TX/RX', 1)
    uhd_usrp_source_0.set_auto_dc_offset(False, 1)
    uhd_usrp_source_0.set_auto_iq_balance(False, 1)
    uhd_usrp_source_0.set_center_freq(center_freq, 2)
    uhd_usrp_source_0.set_gain(rx_gain, 2)
    uhd_usrp_source_0.set_antenna('RX2', 2)
    uhd_usrp_source_0.set_auto_dc_offset(False, 2)
    uhd_usrp_source_0.set_auto_iq_balance(False, 2)

    time.sleep(2)

    print('''\nTo collect a data set press the \'c\' key.
To activate a switch press the \'s\' key.
To end this session press \'q\'.'''
    )
    # while(True)
    user_input = raw_input("Enter your command: ")
    junk = uhd_usrp_source_0.finite_acquisition_v(1000) # first call returns nothing
    data1 = uhd_usrp_source_0.finite_acquisition_v(1000)
    data2 = uhd_usrp_source_0.finite_acquisition_v(1000)
    # print(junk)
    # print(data1)
    pdb.set_trace()
    ##
    # Put switching code here
    #
    ##
    with open('filename1.csv','wb') as out:
        csv_out=csv.writer(out)
        for samp in data1:
            csv_out.writerow([samp.real,samp.imag])

    time.sleep(1)
    data2 = uhd_usrp_source_0.finite_acquisition_v(1)
    print(data2)

    with open('filename2.csv','wb') as out:
        csv_out=csv.writer(out)
        for samp in data2:
            csv_out.writerow([samp.real,samp.imag])

if __name__ == '__main__':
    main()

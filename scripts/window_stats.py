#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
from hirise_assembly import HiriseAssembly
import random
import re
if __name__=="__main__":
     import sys
     import argparse

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug',default=False  ,action="store_true",help="Turn on debugging ouput")
#     parser.add_argument('-L','--layout',default=False ,help="A file containing a layout of contigs.")
     parser.add_argument('-T','--test',  default=False ,help="Compute stats for an explicitly specified region.")
     parser.add_argument('-B','--bin',  default=False ,action="store_true",help="Bin reference sequences.")
     parser.add_argument('-N','--nsamples',default=10000,type=int ,help="Number of winows to sample.")
     parser.add_argument('-w','--window',default=1000, type=int ,help="Size of windows to examine.")
     parser.add_argument('-W','--binwindow',default=100000, type=int ,help="Size of windows to break reference contigs into.")
     parser.add_argument('-i','--infile',default=False ,help="Filename for serialised assembly input file.")
#     parser.add_argument('-o','--outfile',default=False,help="Filename for writing a list of segments on the raw contigs to mask for being promiscuous in linking.")

     args = parser.parse_args()

     if args.infile:
          asf = HiriseAssembly()
          asf.load_assembly(args.infile)
     n=0

     asf.binsize = args.binwindow

     while n<args.nsamples:
          try:
               c,a,b=asf.random_window(wlen=args.window)
               if args.test:
                    m=re.match("^(.*):(\d+)-(\d+)$",args.test)
                    c=m.group(1)
                    a=int(m.group(2))
                    b=int(m.group(3))
               asf.window_stats(c,a,b,debug=args.debug,bins=args.bin)
               n+=1
          except Exception as e:
               print("#",e,c,a,b)
               raise e
          if args.test: break


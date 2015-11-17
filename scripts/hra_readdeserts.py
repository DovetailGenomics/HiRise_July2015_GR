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
#import struct
#import hashlib
import re    

if __name__=="__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d','--debug',default=False,action="store_true",help="Turn on debugging ouput")
    parser.add_argument('-q','--mapq',default=55,  type=float,help="Minimum map quality threshold.")
    parser.add_argument('-i','--infile',default=False,help="Filename for serialised assembly input file.")
    parser.add_argument('-o','--outfile',default=False,help="Filename for bedfile output.")

    args = parser.parse_args()

    if args.infile:
        hra = HiriseAssembly()
        hra.load_assembly(args.infile)

    if args.outfile:
        of=open(args.outfile,"wt")
    else:
        of=sys.stdout

    hra.read_deserts(mapq=55,outfile=of,min_len=1000)

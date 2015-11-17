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

if __name__=="__main__":
     import sys
     import argparse

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug',default=False  ,action="store_true",help="Turn on debugging ouput")
     parser.add_argument('-L','--layout',default=False ,help="A file containing a layout of contigs.")
     parser.add_argument('-i','--infile',default=False ,help="Filename for serialised assembly input file.")
     parser.add_argument('-o','--outfile',default=False,help="Filename for writing a list of segments on the raw contigs to mask for being promiscuous in linking.")

     args = parser.parse_args()

     if args.infile:
          asf = HiriseAssembly()
          asf.load_assembly(args.infile)
     
     if args.outfile:
          f=open(args.outfile,"wt")
          for contig in asf.contigs_iter():
               f.write("{}\t{}\n".format(contig,asf.contig_length(contig)))




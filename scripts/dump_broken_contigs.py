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
from p2fa import fastaWriter

if __name__=="__main__":
     import sys
     import argparse

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug',default=False  ,action="store_true",help="Turn on debugging ouput")
     parser.add_argument('-L','--layout',default=False ,help="A file containing a layout of contigs.")
     parser.add_argument('-i','--infile',default=False ,help="Filename for serialised assembly input file.")
     parser.add_argument('-f','--fasta', default=False ,help="Filename for original contigs fasta.")
     parser.add_argument('-o','--outfile',default=False,help="Filename for writing a list of segments on the raw contigs to mask for being promiscuous in linking.")

     args = parser.parse_args()

     if args.infile:
          asf = HiriseAssembly()
          asf.load_assembly(args.infile)

     asf.ocontig_fasta = args.fasta

     
     
     if args.outfile:
          outfasta=fastaWriter(args.outfile)
          for contig in asf.contigs_iter():
               outfasta.next(contig)
               outfasta.write(asf.get_seq(contig))

          outfasta.flush()
          outfasta.close()



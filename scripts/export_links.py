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
import struct
import hashlib
    

if __name__=="__main__":
     import sys
     import argparse

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug',default=False,action="store_true",help="Turn on debugging ouput")
     parser.add_argument('-C','--nchunks' ,default=1,type=int,help="Number of chunks.")
     
     parser.add_argument('-m','--mask'   ,default=[], action="append",help="Name of a file containing regions of the input assembly to mask out")
     parser.add_argument('-c','--chunk',default=0,  type=int,help="This chunk.")
     parser.add_argument('-K','--contig',default=False,help="Just this contig.")
     parser.add_argument('-q','--mapq',default=10,  type=float,help="Minimum map quality threshold.")
     parser.add_argument('-i','--infile',default=False,help="Filename for serialised assembly input file.")
     parser.add_argument('-o','--outfile',default=False,help="Filename for writing a list of segments on the raw contigs to mask for being promiscuous in linking.")
# -m 2 -w 1000 -M $( cat {input.threshold} ) 
     args = parser.parse_args()

     if args.infile:
          asf = HiriseAssembly()
          asf.load_assembly(args.infile)
     
     for segments_file in args.mask:
         asf.add_mask_regions(filename=segments_file)
         asf.merge_masked_regions()

     
     if args.contig:
          contig_iter=iter([args.contig])
     else:
          contig_iter = asf.ocontigs_iter()

     if args.outfile:
          of=open(args.outfile,"wt")
     else:
          of=sys.stdout
     for ocontig in contig_iter:
          chunk = struct.unpack("<L", hashlib.md5(ocontig.encode("utf-8")).digest()[:4])[0]%args.nchunks
          if (not args.contig) and (not chunk == args.chunk): continue
          links={}
          asf.get_links([ocontig],skipI=True,mapq=args.mapq,links=links,contigs=False,raw=False,debug=args.debug)
          for c1,c2 in links.keys():
               l1,l2=asf.contig_length(c1),asf.contig_length(c2)
               l=links[c1,c2]
               n=len(l)
               of.write("\t".join(map(str,[c1,c2,l1,l2,n,l]))+"\n")
     of.close()

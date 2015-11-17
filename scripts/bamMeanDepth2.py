#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
import sys
import argparse
import pysam

if __name__=="__main__":

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug',default=False,action="store_true",help="Turn on debugging ouput")
     parser.add_argument('-p','--progress',default=False,action="store_true",help="Print progress info")
     parser.add_argument('-b','--bamfile',default=False,help="Bam file.")
     parser.add_argument('-B','--bedfile',default=False,help="Bed file.")
     parser.add_argument('-q','--mapq',default=10,type=int,help="Mapq.")
     parser.add_argument('-o','--outhist')

     args = parser.parse_args()
     if args.progress: print("#",args)

     bamfile=pysam.Samfile(args.bamfile)

     h=bamfile.header
     seqs=h['SQ']

     slen=[ s['LN'] for s in seqs ]
     snam=[ s['SN'] for s in seqs ]

     llen={}
     for i in range(len(snam)):
          llen[snam[i]]=slen[i]
     
     out=0
     in1=1
     in2=2

     state=out

     h={}

     last_scaffold=False
     last_x="-"
     start_x=0

#     for c in mycontigs:
#          if llen[c]<1000: continue
#          db=[]
#          nr=0

     tl=0
     nr=0

     for aln in bamfile.fetch(until_eof=True):
          if aln.is_duplicate     : continue
          if aln.mapq < args.mapq : continue
          scaffold=bamfile.getrname(aln.tid)
          c=scaffold
          x = aln.pos
          if x < 500 or x > llen[c]-500: continue
          nr+= 1
          tl+= aln.query_alignment_length
#          print("#",aln)

          if last_scaffold and not scaffold==last_scaffold:
               if not llen[last_scaffold]<1000: 
                    c=last_scaffold
                    print(c,tl/(llen[c]-999),llen[c],nr,sep="\t")
#               else:
#                    print("#",last_scaffold)
               nr=0
               tl=0
          last_scaffold=scaffold

     if not llen[last_scaffold]<1000: 
          c=last_scaffold
          print(c,tl/(llen[c]-999),llen[c],nr,sep="\t")


               




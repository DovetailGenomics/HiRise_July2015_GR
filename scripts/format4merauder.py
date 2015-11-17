#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import map
from builtins import range
from past.utils import old_div
import sys
import argparse


if __name__=="__main__":

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug'   ,default=False,action="store_true",help="Turn on debugging ouput")
     parser.add_argument('-k','--kmersize',default=51,type=int,help="Kmer size -- gaps with shorter kmers a skipped")
     
#     parser.add_argument('-p','--progress',default=False,action="store_true",help="Print progress info")
#     parser.add_argument('-L','--length',default=False,type=int,help="File containing lenghts.")

     args = parser.parse_args()
#     if args.progress: print("#",args)

     data=[]
     lastp=(0,0)
     lasts=-1
     lastk1=0
     lastk2=0
     while True:
          l=sys.stdin.readline()
          if not l: break
          if l[0]=="#": continue
          c=l.strip().split()
          if not len(c)>=8: 
               print("##",c)
               exit(0)
          k1=c[7]
          k2=c[8]
          p=(c[3],c[4])
          if not p==lastp:
               if not lastp==(0,0):
                    if len(lastk1)==args.kmersize and len(lastk2)==args.kmersize:
                         print("\t".join(map(str,[ lasts, lastp[0], lastk1, lastp[1], lastk2, 5,200, "\t".join(data) ])))
               data=[]

          data.append(c[11])

          lasts=c[1]
          lastk1=k1
          lastk2=k2
          lastp=p

     if not lastp==(0,0):
          if len(lastk1)==args.kmersize and len(lastk2)==args.kmersize:
               print("\t".join(map(str,[ lasts, lastp[0], lastk1, lastp[1], lastk2, 5,200, "\t".join(data) ])))

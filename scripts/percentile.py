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
import gc

if __name__=="__main__":

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug',default=False,action="store_true",help="Turn on debugging ouput")
     parser.add_argument('-p','--percentile',default=99.0 ,type=float, help="Percentile")
     parser.add_argument('-m','--minvalue',default=False ,type=float, help="Don't return a value lower than X")
     parser.add_argument('-L','--length'    ,default=False,type=int   ,help="File containing lenghts.")

     args = parser.parse_args()

     buffer=[]

     gc.disable()
     while True:
          l=sys.stdin.readline()
          if not l: break
          if l[0]=="#": continue
          x=float(l.strip().split()[0])
          buffer.append(x)
     gc.enable()
     
     buffer.sort()
     l=len(buffer)
     i=int(l*args.percentile/100)
     if args.minvalue:
          print(max(int(buffer[i]),args.minvalue))
     
     else:
          print(int(buffer[i]))

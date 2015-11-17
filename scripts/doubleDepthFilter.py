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

     args = parser.parse_args()

     b=[]
     total_l=0
     for l in sys.stdin:
          c=l.strip().split()
          if len(c)<2: break
          b.append( (float(c[1]),c[0],int(c[2])) )
          total_l += int(c[2])
     b.sort()
     medi=0
     rs=0
     while rs < total_l/2:
          rs+=b[medi][2]
          medi+=1
#     medi=b[int(len(b)/2)]
     #     print("#",medi)

     threshold_depth = b[medi][0]*1.5

     for x,s,l in b:
          if x>threshold_depth:
#               print(s,x,sep="\t")
               print(s)



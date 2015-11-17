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


if __name__=="__main__":

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug',default=False,action="store_true",help="Turn on debugging ouput")
     parser.add_argument('-p','--progress',default=False,action="store_true",help="Print progress info")

     args = parser.parse_args()
     if args.progress: print("#",args)

     x=0
     state="OUT"
     while True:
          l=sys.stdin.readline()
          if not l: break
          if l.startswith('>'): 
               name = l[1:].strip().split()[0]
#               print(">>",name)
               x=0
               state="OUT"
          else:
               s=l.strip()

               i=x
#               print(name,i,s)
               while i < x+len(s):
                    
                    if s[i-x] in ["N","n"] and state=="OUT":
#                         print(s[i-x])
                         state="IN"
                         xstart=i
                    if ( not s[i-x] in ["N","n"]) and state=="IN":
                         state="OUT"
                         print("GAP",name,xstart,i,i-xstart,"gap",sep="\t")
#                         xstart=i
                         

                    i+=1

               x+=len(s)
     if state=="IN":
          print("GAP",name,xstart,x,x-xstart,"gap",sep="\t")
          

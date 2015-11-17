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
     parser.add_argument('-b','--breaks',default=False ,help="File containing breaks")
     parser.add_argument('-i','--infile',default=False ,help="Input layout in p: lines format.")
     parser.add_argument('-o','--outfile',default=False,help="Filename for output.")

     args = parser.parse_args()


     if args.infile:
          asf = HiriseAssembly()
          asf.load_assembly(args.infile)


#     asf = HiriseAssembly()
#     asf.load_playout(args.infile)

     breaks=[]
     scores={}
     for l in open(args.breaks):
          if l[0]=="#": continue
#Scaffold102239 741 1379 5097 -5.421961663655971
          scaffold,a,b,slen,score = l.strip().split()
          a=int(a) #start
          b=int(b) #end
#          c=int(c) #lowpoint
          score=float(score)
          c=int((a+b)/2)
          breaks.append((scaffold,a,b,c))
          scores[scaffold,a,b,c]=score

#     asf.add_breakpoints_ranges(breaks,debug=args.debug,scores=scores)
     asf.add_breakpoints_ranges(breaks,debug=args.debug)
     asf.validate()

     if args.outfile: asf.save_assembly( args.outfile )

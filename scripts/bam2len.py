#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from past.utils import old_div
import pysam
import sys

def log(s):
    sys.stderr.write(  "%s\n" %(s) )

if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--bamfile',required=True)
    parser.add_argument('-d','--debug',default=False,action="store_true")
    parser.add_argument('-p','--progress',default=False,action="store_true")
    parser.add_argument('-c','--cutoff',default=(old_div(200.0,1.0e13)),type=float)

    nodes={}

    args = parser.parse_args()
    if args.debug:
        args.progress=True

    if args.progress: log( str(args) )
  
    sam=pysam.Samfile(args.bamfile)

    if args.progress: log( "opened %s" % args.bamfile )


    h=sam.header
    seqs=h['SQ']

    if args.progress: log(  "%s %d"% ( "number of sequences in the reference:", len(seqs) )   )

    slen=[ s['LN'] for s in seqs ]
    snam=[ s['SN'] for s in seqs ]

    if args.progress: log(  "built length and name map arrays" )   


    for i in range(len(seqs)):
        print(snam[i],slen[i])

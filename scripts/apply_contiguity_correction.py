#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import str
from past.utils import old_div
import sys
#import networkx as nx

if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-c','--correction',required=True)
#    parser.add_argument('-G','--genomesize',type=float)
    parser.add_argument('-m','--max',type=int,default=1000000)
    parser.add_argument('-d','--debug',action="store_true")
    parser.add_argument('-p','--progress',action="store_true")

    args = parser.parse_args()
    if args.debug:
        args.progress=True

    if args.progress: print("#",str(args)) 
    sys.stdout.flush()
    correction={}
    G=0.0
    if args.correction:
        f = open(args.correction)
        while True:
            l = f.readline()
            if not l: break
            if l[0]=="#": continue

            c=l.strip().split()
            correction[int(c[0])]=float(c[1])

        f.close()

    while True:
        l=sys.stdin.readline()
        if not l: break
        c=l.strip().split()
        print(c[0],old_div(float(c[1]),correction.get(int(float(c[0])),1.0)))


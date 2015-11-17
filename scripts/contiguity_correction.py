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
import sys
#import networkx as nx

if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-l','--lengths')
    parser.add_argument('-G','--genomesize',type=float)
    parser.add_argument('-m','--max',type=int,default=1000000)
    parser.add_argument('-d','--debug',action="store_true")
    parser.add_argument('-p','--progress',action="store_true")

    args = parser.parse_args()
    if args.debug:
        args.progress=True

    if args.progress: print("#",str(args)) 
    sys.stdout.flush()
    ll=[]
    G=0.0
    if args.lengths:
        f = open(args.lengths)
        print("#load",args.lengths)
        while True:
            l = f.readline()
            if not l: break
            if l[0]=="#": continue

            c=l.strip().split()
            ll.append(int(c[1]))
            G+=int(c[1])
        f.close()

    print("#done loading lengths.")
    ll.sort()

    if args.genomesize:
        G=args.genomesize

    lc=[]
    rs=0.0
#    G=sum(ll)

    for l in ll:
        rs+=l
        lc.append(rs)

    nc=len(ll)
    import bisect

    if args.debug:
        for i in range(len(ll)):
            print(i,ll[i],lc[i])

    for i in range(args.max):


        j=bisect.bisect_left( ll , i )
        if j>=len(lc): break
        #for l in ll:
        #    if l>i:
        #        x+=l-i
        #    else:
        #        break
        print(i,old_div((G+(nc-j)*i-lc[j]),G)) #,j,ll[j],G,nc-j, lc[j]>i



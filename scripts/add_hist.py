#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
from builtins import range
import sys

if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-m','--max',type=int,default=1000000)
    parser.add_argument('-d','--debug',action="store_true")

    parser.add_argument('files', metavar='FILENAME', nargs='+',
                        help='histogram files')

    args = parser.parse_args()
    if args.debug:
        args.progress=True

    cs=[]

    s1={}
    s2={}
    
    min_i=1e6
    max_i=0

    for fn in args.files:
#        c1={}
        f = open(fn)
        while True:
            l = f.readline()
            if not l: break
            if l[0]=="#": continue

            c=l.strip().split()
 #           c1[int(c[0])]=(int(c[1]),int(c[2]))
            x=int(c[0])
            if x>args.max: break
            a=int(c[1])
            s1[x]=s1.get(x,0)+a
            if len(c)>2:
                b=int(c[2])
                s2[x]=s2.get(x,0)+b
            if x>max_i: max_i=x
            if x<min_i: min_i=x
        f.close()

    for i in range(min_i,max_i+1):
        print(i, s1.get(i,0), s2.get(i,0))

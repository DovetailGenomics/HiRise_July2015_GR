#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#


## Double threshold

from __future__ import print_function
import sys
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('-t', "--thresholds", required=True)
parser.add_argument('-d','--debug',action="store_true")

args = parser.parse_args()

with open(args.thresholds, "rt") as thr_handle:
    t1, t2 = [float(l.rstrip()) for l in thr_handle]

out=0
in1=1
in2=2

state=out


last_scaffold="-"
last_x="-"
start_x=0
while True:
    l= sys.stdin.readline()
    if not l: break
    c = l.strip().split()
    scaffold,x,y = c[0],int(c[1]),int(c[2])
    if args.debug: print("#",scaffold,x,y,state,start_x)

    if scaffold==last_scaffold:
        if state == in2 and y < t1:
            print(scaffold,start_x,x,"deep")
            state=out
            continue
    else: 
        if state==in2:
            print(last_scaffold,start_x,last_x,"deep") #, "#toend")
        state=out
        last_scaffold=scaffold
        last_x=x
        continue

    if state==out and y >= t1: 
        state=in1
        start_x=x
        
    if y >= t2: 
        state=in2
    

    if y< t1: state=out
    last_scaffold=scaffold
    last_x=x


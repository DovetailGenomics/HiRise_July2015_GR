#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#


from __future__ import print_function
import sys

def process_buffer(b):
    b.sort(key=lambda x: (x[0],-x[1],x[2]))
    rs=0
    state=0
    for bb in b: 
        rs+=bb[1]
        x=bb[0]
        #print bb,rs
        if state==0 and rs>0:
            begin=x
            state=1
        elif state==1 and rs==0:
            state=0
            print(bb[2],begin,x-1,bb[3],bb[4])
#        lastx=x

if __name__=="__main__":

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-b','--buffer',default=1000,type=int)
    parser.add_argument('-m','--minlen',default=2000,type=int)
    parser.add_argument('-t','--thresh',default=1.0,type=float)
    parser.add_argument('-n','--dryrun',default=False,action="store_true")

    parser.add_argument('-w','--window',   default=50   ,type=int)
    parser.add_argument('-d','--debug',    default=False,action="store_true")

    parser.add_argument('-p','--progress', default=False,action="store_true")

    args = parser.parse_args()

#    print args

    b=[]
    lastc=False
    while True:
        l=sys.stdin.readline()
        if not l: break
        if l[0]=="#": continue
#        print l

        c=l.strip().split()
        scaffold,startx,endx,length,minl = c[0],int(c[1])+1,int(c[2])+1,int(c[3]),float(c[4])
        if length        < args.minlen: continue
        if startx        < args.buffer: continue
        if (length-endx) < args.buffer: continue
        if minl > args.thresh         : continue

        if lastc and not c[0]==lastc:
            process_buffer(b)
            b=[]


        sl=int(c[3])
        b.append((max(0   ,int(c[1])-args.window), 1,c[0],sl,float(c[4])))
        b.append((min(sl-1,int(c[2])+args.window),-1,c[0],sl,float(c[4])))

        lastc=c[0]
    process_buffer(b)




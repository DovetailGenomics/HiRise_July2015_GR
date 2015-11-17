#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
from builtins import map
import sys
import argparse

#Scaffold1       CONTIG1 +isotig666668   1       581   

parser = argparse.ArgumentParser()

#parser.add_argument('-i','--input')
parser.add_argument('-d','--debug',default=False,action="store_true")
parser.add_argument('-p','--progress',default=False,action="store_true")
parser.add_argument('-L','--length',default=False,type=int)

args = parser.parse_args()
if args.progress: print("#",args)

scaffolds={}
last_contig=-1
contig_strand={}
contigx={}
ocontig_contigs={}
contig_length={}
last_scaffold=-1
contig_scaffold={}
lastx=-1
#scaffold_gaps={}

while True:
        l=sys.stdin.readline()
        if not l: break
        c=l.strip().split()
        if not c[0]=="p:": continue


        scaffold,end = c[1],c[2]
        x=int(float(c[3]))
        scaffolds[scaffold]=1
        contig = end[:-2]
        name_pieces=contig.split("_")
        base = name_pieces[-1]
        ocontig = "_".join(name_pieces[:-1])

        if not last_contig==contig:
             contig_strand[contig]=1 if (end[-1]=="5") else -1
#        print contig_strand[contig],contig
        if end[-1]=="5": 
            contigx[contig]=x
            ocontig_contigs[ocontig] = ocontig_contigs.get(ocontig,[])+[(int(base)-1,contig)]
            ocontig_contigs[ocontig].sort()

        length=int(c[7])
        contig_length[contig]=length
        blast_coord = int(c[5])

        if not last_scaffold==scaffold:
             contign=1

        if last_scaffold==scaffold and not last_contig==contig: 
             pass
             print("\t".join(map(str,[scaffold,"GAP",x-lastx-1,50])))
#             print "gap",lastx,x,scaffold
#            scaffold_gaps[scaffold]=scaffold_gaps.get(scaffold,[]) + [(lastx,x)]
#            gap2ends[(scaffold,lastx,x)] = (last_end,end)
        elif last_scaffold==scaffold: # not last_scaffold==scaffold:
             st="-" if contig_strand[contig]==-1 else "+"
             print("\t".join(map(str,[scaffold,"CONTIG{}".format(contign),"{}{}".format(st,contig),lastx,x,20.0])))
             contign+=1
             pass
            #            print "gaps:",scaffold,scaffold_gaps.get(last_scaffold,[])

        blast_chr = c[4]

        contig_scaffold[contig] =scaffold
        contig_scaffold[ocontig]=scaffold

        bh=eval(" ".join(c[8:]))
        if bh:
            unc =abs(length- abs(int(bh[3])-int(bh[4])))

            if bh[2]=="+" and end[-1]=="5": blast_coord -= int(bh[5])-1
            if bh[2]=="+" and end[-1]=="3": blast_coord += length-int(bh[6])-1
            if bh[2]=="-" and end[-1]=="5": blast_coord += int(bh[5])-1
            if bh[2]=="-" and end[-1]=="3": blast_coord -= length-int(bh[6])-1

        else:
            unc=0

        last_scaffold=scaffold
        last_contig=contig
        lastx=x
        last_end=end


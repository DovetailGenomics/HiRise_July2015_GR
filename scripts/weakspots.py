#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
from builtins import range
import sys
import mapper
from collections import deque
#import heapq
from itertools import chain
import bisect
import re
import chicago_edge_scores as ces
import math
from bamtags import BamTags
import pysam

def log(s):
    sys.stderr.write(  "%s\n" %(s) )

if __name__=="__main__":
    import sys

#    print " ".join(sys.argv)

    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c','--chunk',required=True)
    parser.add_argument('-b','--bamfiles',required=True)
    parser.add_argument('-i','--ignoreDuplicates'   ,default=True,action="store_false")
    parser.add_argument('-q','--mapq',default=10    ,type=int,help=" ")
    parser.add_argument('-t','--trim',default=100   ,type=int,help=" ")
    parser.add_argument('-L','--minlen',default=1000,type=int,help=" ")

    parser.add_argument('-A','--end_threshold',  default=25.0     , type=float, help=" ")
    parser.add_argument('-B','--support_cutoff', default=10.0     , type=float, help=" ")
    parser.add_argument('-R','--range_cutoff',   default=200000.0 , type=float, help=" ")

    parser.add_argument('-d','--debug',default=False,action="store_true")
    parser.add_argument('-p','--progress',default=False,action="store_true")
    parser.add_argument('-M','--set_insert_size_dist_fit_params')

    args = parser.parse_args()
    nodes={}

    range_cutoff=args.range_cutoff
    t1=args.end_threshold
    t2=args.support_cutoff

#    recombinings2.max_n_breaks=args.max_n_breaks

    fmodel=open( args.set_insert_size_dist_fit_params )
    contents = fmodel.read()
    try:
        fit_params=eval(contents)
    except:
        "couldn't deal with option", args.param
    fmodel.close
    ces.set_exp_insert_size_dist_fit_params(fit_params)

    cache_cutoff = 150000
    log_score_cache={}

    bamfiles = args.bamfiles.split(",")

    contigs = []
    chunkf=open(args.chunk,"rt")
    while True:
        l=chunkf.readline()
        if not l: break
        contigs.append(l.strip())
    chunkf.close()

    oname={}
    bams=[]
    for bamfile in bamfiles:
        bams.append( pysam.Samfile(bamfile,"rb") )
        oname[bams[-1]]=bamfile
    h=bams[0].header
    seqs=h['SQ']

    slen=[ s['LN'] for s in seqs ]
    snam=[ s['SN'] for s in seqs ]

    llen={}
    for i in range(len(snam)):
        llen[snam[i]]=slen[i]

    buff=0 #args.trim
    for contig in contigs:
        if args.debug: print("#{}".format(contig))
        if (not contig in llen) or  llen[contig]<args.minlen: continue
        spans=[]
        edges=[]
        nb=[]
        lll = llen[contig]
#        print "#",contig
        for sam in bams:
            for aln in sam.fetch(region=contig):
                if (aln.tid==aln.rnext) and (aln.pos<aln.pnext) and ((not args.ignoreDuplicates) or (not aln.is_duplicate)) and (aln.mapq >= args.mapq) and (BamTags.mate_mapq(aln) >= args.mapq):
#                   spans.append( tuple( [aln.pos, aln.pos+aln.tlen,aln.mapq,oname[sam]] )  )
#                   if aln.tlen > 2*buff and aln.tlen < range_cutoff:
                   if aln.tlen > 2*buff: # and aln.tlen < range_cutoff:
                       ll=ces.model.lnF(aln.tlen)
#                       ll=get_score( aln.tlen )                   
                       edges.append( tuple([aln.pos + buff         , ll]) )
                       edges.append( tuple([min(lll,aln.pos+aln.tlen-buff)  ,-ll]) )

                       nb.append( tuple([aln.pos + buff         , 1]) )
                       nb.append( tuple([min(lll,aln.pos+aln.tlen-buff)  ,-1]) )

        edges.sort()
        nb.sort()
        rs=0
        n=0
        tripped=False
        last=0
        state=0
        stretch_start=0
        low_point=0.0
        for i in range(len(edges)):
            rs+=edges[i][1]
            n+=nb[i][1]
            x=edges[i][0]

            try:
                score=ces.model.cutScore(llen[contig],x,n,rs)
            except Exception as e:
                for ii in range(len(edges)):
                    print(ii,edges[ii])
                print("exception at i={}, x={}, len={}".format(i,x,llen[contig]))
                raise e
            if score>t1:
                tripped=True
            if tripped and score<t2 and state==0:
                stretch_start=x
                state=1
                low_point =score
            if state==1 and score>t2:
                print(contig,stretch_start,x,llen[contig],low_point,"break")
                state=0
            if state==1:
                if score<low_point:
                    low_point=score
            if args.debug: print("dd:",contig,edges[i][0],rs,state,x,stretch_start,score,n)

            last=edges[i][0]

    exit(0)



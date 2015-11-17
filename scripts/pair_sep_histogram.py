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
import networkx as nx
import sys
from bamtags import BamTags


def log(s):
    sys.stderr.write(  "%s\n" %(s) )

if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--bamfile',required=True)
    parser.add_argument('-L','--contiglist')
#    parser.add_argument('-s','--savegraph',default=False)
#    parser.add_argument('-S','--restoregraph',default=False)
    parser.add_argument('-H','--head',type=int,default=False)
    parser.add_argument('--insertHist',default=False,action="store_true")
    parser.add_argument('-d','--debug',default=False,action="store_true")
    parser.add_argument('-m','--minlength',default=1000,type=int)
    parser.add_argument('-q','--min_qual',default=10,type=int)
#    parser.add_argument('-e','--endsonly',default=False,action="store_true")
#    parser.add_argument('-E','--endwindow',default=50000,type=int)
#    parser.add_argument('-D','--dotfile',default=False)
    parser.add_argument('-p','--progress',default=False,action="store_true")
    parser.add_argument('-c','--cutoff',default=(old_div(200.0,1.0e13)),type=float)
    parser.add_argument('-j', '--junctions', default=False, action="store_true", help="Only count reads with junctions")

    g = nx.Graph()
    g2 = nx.Graph()

    nodes={}

    args = parser.parse_args()
    if args.debug:
        args.progress=True

    contigs=[]
    if args.contiglist:
        f=open(args.contiglist)
        while True:
            l=f.readline()
            if not l: break
            c=l.strip().split()
            if l[0]=="#": continue
            contigs.append(c[0])

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
        if slen[i]>args.minlength:
            nodes[i]=1

    n=0
    n2=0
    ne=0
    hist ={}
    hist2={}
    nr=0
    if args.progress: log("about to iterate over bamfile alignments")

    if args.contiglist:

        def alngenerator():
            for c in contigs:
                sys.stderr.write("{}\n".format(c))
                for aln in sam.fetch(region=c):
                    yield aln

        alniter = alngenerator()
    else:
        alniter = sam.fetch(until_eof=True)

    for aln in alniter:
#        print aln.pos, dir(aln)
        nr+=1
        if nr%200000 == 0:
            sys.stderr.write("%d\n"%nr)

        i,j = aln.tid, aln.rnext
        if args.debug: print([i,j,int(aln.mapq),aln.tlen,aln.is_duplicate])

        if int(aln.mapq) < args.min_qual: continue
        if BamTags.mate_mapq(aln) < args.min_qual: continue
        if aln.is_duplicate : continue
        if args.junctions and BamTags.junction(aln) != "T": continue

        if i==j:
            if (slen[i]<args.minlength) : continue
            if aln.tlen < 0: continue
            hist[aln.tlen] = hist.get(aln.tlen,0)+1
            n+=1
        else:
            if (slen[i]<args.minlength) or (slen[j]<args.minlength) : continue
            if not aln.is_read1: continue
            x = min( aln.pos , slen[i] - aln.pos ) + min( aln.pnext,slen[j]-aln.pnext)
            hist2[x] = hist2.get(x,0)+1 
            n2+=1

        if args.head and (n+n2)>args.head:
            break

#    k = list(hist.keys())
#    k.sort()

    nnn = max( list(hist.keys())+ list(hist2.keys()) )
    for i in range(nnn+1):
        print(i,hist.get(i,0),hist2.get(i,0))

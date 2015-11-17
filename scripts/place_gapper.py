#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
from builtins import map
from builtins import range

import sys

import pysam

import mapper
from collections import deque

from itertools import chain
import bisect
import re

import math
import hashlib
import struct
from mapperlite import MapperLite


from bamtags import BamTags

tr = str.maketrans("ACTGactg","TGACtgac")

def rc(s):
    r=s[::-1]
    rrc = r.translate(tr)
    return(rrc)

def log(s):
    sys.stderr.write(  "%s\n" %(s) )

if __name__=="__main__":
    import sys

#    print " ".join(sys.argv)

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--my_chunk',type=int,default=1)
    parser.add_argument('-C','--nchunks',type=int,default=32)
    parser.add_argument('-b','--bamfiles',required=True)
#    parser.add_argument('-i','--ignoreDuplicates',default=True,action="store_false")
    parser.add_argument('-q','--mapq',default=10,type=int)
    parser.add_argument('-t','--trim',default=100,type=int)
    parser.add_argument('-L','--minlen',default=1000,type=int)
    parser.add_argument('-K','--kmer',default=51,type=int)

    parser.add_argument('-f','--fasta',required=True)
    parser.add_argument('-r','--region',required=False,default=False)
    parser.add_argument('-A','--end_threshold',default=25.0,type=float)
    parser.add_argument('-B','--support_cutoff',default=10.0,type=float)
    parser.add_argument('-R','--range_cutoff',default=200000.0,type=float)
    parser.add_argument('-n','--pn',type=float)

    parser.add_argument('--direct',default=False,action="store_true")
    parser.add_argument('--twopass',default=False,action="store_true")
    parser.add_argument('-d','--debug',default=False,action="store_true")
    parser.add_argument('-p','--progress',default=False,action="store_true")

    args = parser.parse_args()
    nodes={}

    range_cutoff=args.range_cutoff
    t1=args.end_threshold
    t2=args.support_cutoff

#    recombinings2.max_n_breaks=args.max_n_breaks

    bamfiles = args.bamfiles.split(",")

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

#
# Read in a hirise layout (from "^p:" lines)
#
    mapper = MapperLite()
    mapper.load_layout(sys.stdin)

    my_scaffolds={}
    scaffold_hashes={}
    for s in list(mapper.scaffolds.keys()):
        scaffold_hashes[s]=struct.unpack("<L", hashlib.md5(s.encode("utf-8")).digest()[:4])[0]%args.nchunks
        if scaffold_hashes[s]==args.my_chunk:
            my_scaffolds[s]=1
            if args.debug: print("my scaffold:",s)

    sin=False
    contig_seq={}
    fh = open(args.fasta)
    name=""
    while True:
        l=fh.readline()
        if not l: break
        if l[0]==">":
            c=l.strip().split()
            name= c[0][1:]
            sc = mapper.contig_scaffold.get(name)
            if args.debug: print("#",name,sc, sc in my_scaffolds,type(sc), my_scaffolds.get(sc),scaffold_hashes.get(sc))
            if mapper.contig_scaffold.get(name) in my_scaffolds:
                contig_seq[name]=""
                sin=True
            else:
                sin=False
        else:
            if sin:
                contig_seq[name]+=l.strip()
    insert_size=500.0


    second_pass_reads = {}

    import struct
    nh=0
    for b in bams:
        for aln in b.fetch(until_eof=True,region=args.region):
            if aln.is_duplicate: continue

            contig = snam[aln.tid]
            ncontig= snam[aln.rnext]

            scaffold ,z1a,z2a,z3a,cdz= mapper.mapCoord( contig, aln.pos,     aln.pos+1 ) 
            nscaffold,z2p,x2p,z3p,cdz=mapper.mapCoord(ncontig, aln.pnext,aln.pnext+1 ) 
 
            if not ((scaffold in my_scaffolds) or ((not args.twopass) and nscaffold in my_scaffolds)): continue

            if args.debug: print("xx",contig,mapper.ocontig_contigs.get(contig), ncontig,mapper.ocontig_contigs.get(ncontig))

            if aln.mapq>=10.0 and scaffold in my_scaffolds : # we might want this read, if it hangs into a gap, or its pair in the next pass

                if args.twopass and not args.direct :  # see if we want to flag this read's pair for extraction later, in the second pass.  if so, store in second_pass_reads
                    template_id = aln.query_name
                    sc,xx,yy,st = scaffold,z1a,z2a,z3a
                    
                    scaffold_strand = st if not aln.is_read2 else -1*st
                    if scaffold_strand == -1:
                        xx -= insert_size
                        yy -= insert_size
                    else:
                        xx += insert_size
                        yy += insert_size
                    gaps = mapper.hits_gaps(sc,xx-400,yy+400,args.debug)
                    if args.debug: print("#gaps:",sc,xx,yy,len(gaps),gaps)
                    if gaps:
#                        seq = aln.seq
#                        qal = aln.qual
#                        if scaffold_strand==1:
#                            seq = rc(seq)
#                            qal = qal[::-1]
                        if args.debug: print("\t".join(map(str,["gaps: p",len(gaps),ncontig,template_id,aln.is_read2,aln.pnext,sc,xx,yy,st,gaps,[mapper.gap2ends[sc,xxx,yyy] for xxx,yyy in gaps]])))

                        second_pass_reads[template_id,aln.is_read2]=[scaffold_strand,[]]
                        for xxx,yyy in gaps:
                            c1,c2= mapper.gap2ends[sc,xxx,yyy]
                            if c1[-1]=="3":
                                k1 =    contig_seq[c1[:-2]][-args.kmer:]
                            else:
                                k1 = rc(contig_seq[c1[:-2]][:args.kmer ])
                            if c2[-1]=="5":
                                k2 =    contig_seq[c2[:-2]][:args.kmer]
                            else:
                                k2 = rc(contig_seq[c2[:-2]][-args.kmer:])
                            #print("\t".join(map(str,["gP:",sc,int(xx),c1,c2,c1[-1],c2[-1],k1,k2,k1 in seq,k2 in seq,"{}:{}".format(seq,qal)])))
                            second_pass_reads[template_id,aln.is_read2][1].append( ["gP:",sc,int(xx),c1,c2,c1[-1],c2[-1],k1,k2,st,aln.is_reverse] )


                sc,xx,yy,st,cdz = mapper.mapCoord( contig, aln.pos,aln.aend )
                if args.debug: print(sc,xx,yy,st)
                gaps = mapper.hits_gaps(sc,xx-50,yy+50)
                if args.debug: print("#",gaps)
                if gaps:
                    seq = aln.seq
                    qal = aln.qual
                    if (st==-1): 
#                    if (st==-1 and not aln.is_reverse) or (st==1 and aln.is_reverse): 
                        seq = rc(seq)
                        qal = qal[::-1]
                    if args.debug: print("\t".join(map(str,["gaps: s",contig,aln.pos,sc,xx,yy,st,gaps,[mapper.gap2ends[sc,xxx,yyy] for xxx,yyy in gaps]])))
                    for xxx,yyy in gaps:
                        c1,c2= mapper.gap2ends[sc,xxx,yyy]
                        if not c1[:-2] in contig_seq:  print("wft? c1",c1,sc,xx,yy,st,scaffold,nscaffold,scaffold in my_scaffolds,mapper.contig_scaffold[c1[:-2]],type(sc),type(scaffold),scaffold_hashes[scaffold])
                        if not c2[:-2] in contig_seq:  print("wft? c2",c2,sc,xx,yy,st,scaffold,nscaffold,scaffold in my_scaffolds,mapper.contig_scaffold[c2[:-2]],type(sc),type(scaffold),scaffold_hashes[nscaffold])

                        #
                        #  This is not right:  assumes that the two reads hit the contigs on either side of the gap... 
                        #
                        if c1[-1]=="3":
                            k1 = contig_seq[c1[:-2]][-args.kmer:]
                        else:
                            k1 = rc(contig_seq[c1[:-2]][:args.kmer])
                        if c2[-1]=="5":
                            k2 = contig_seq[c2[:-2]][:args.kmer]
                        else:
                            k2 = rc(contig_seq[c2[:-2]][-args.kmer:])
                        print("\t".join(map(str,["gS:",sc,int(xx),c1,c2,c1[-1],c2[-1],k1,k2,k1 in seq, k2 in seq,"{}:{}".format(seq,qal),st,aln.is_reverse])))
                        #if overlaps_gap(aln2scaffoldx(aln)):
                        #pass
                        # this read overlapps a gap?
                        #OPEN:Scaffold278839     isotig431029.5  TTTTTCTCATCAGTCTCCTCATCACAAGTTTTCTATAAGTCCCCA   isotig1560433.5 TTAAATCATAGCAGTACTCATGAAGAGGGTAGTAAAGCACTGGAA   57      14      :       AG

            if (not args.direct) and (not args.twopass) and (BamTags.mate_mapq(aln)>=10.0 and nscaffold in my_scaffolds) : # we might want this read because it's sister maps to one of our scaffolds
                sc,xx,yy,st,cdz = mapper.mapCoord( ncontig,aln.pnext,aln.pnext+1 )
                mate_scaffold_strand = st if not aln.mate_is_reverse else -1*st
                if mate_scaffold_strand == -1:
                    xx -= insert_size
                    yy -= insert_size
                else:
                    xx += insert_size
                    yy += insert_size
                gaps = mapper.hits_gaps(sc,xx-400,yy+400)
                if args.debug: print("#",gaps)
                if gaps:
                    seq = aln.seq
                    qal = aln.qual
                    if mate_scaffold_strand==1:
                        seq = rc(seq)
                        qal = qal[::-1]
                    if args.debug: print("\t".join(map(str,["gaps: p",ncontig,aln.pnext,sc,xx,yy,st,gaps,[gap2ends[sc,xxx,yyy] for xxx,yyy in gaps],seq,qal])))
                    for xxx,yyy in gaps:
                        c1,c2= mapper.gap2ends[sc,xxx,yyy]
                        if c1[-1]=="3":
                            k1 =    contig_seq[c1[:-2]][-args.kmer:]
                        else:
                            k1 = rc(contig_seq[c1[:-2]][:args.kmer ])
                        if c2[-1]=="5":
                            k2 =    contig_seq[c2[:-2]][:args.kmer]
                        else:
                            k2 = rc(contig_seq[c2[:-2]][-args.kmer:])
                        print("\t".join(map(str,["gP:",sc,int(xx),c1,c2,c1[-1],c2[-1],k1,k2,k1 in seq,k2 in seq,"{}:{}".format(seq,qal)])))



                nh+=1
        b.close()

    if args.twopass and not args.direct:    
        bams=[]
        for bamfile in bamfiles:
            bams.append( pysam.Samfile(bamfile,"rb") )
            oname[bams[-1]]=bamfile
        for b in bams:
            for aln in b.fetch(until_eof=True,region=args.region):
                if aln.is_duplicate: continue
                template_id = aln.query_name
                if (template_id, not aln.is_read2) in second_pass_reads:
                    scaffold_strand = second_pass_reads[template_id,not aln.is_read2][0]
                    seq = aln.seq
                    qal = aln.qual
                    if scaffold_strand==1:
                        seq = rc(seq)
                        qal = qal[::-1]

                    for row in  second_pass_reads[template_id,not aln.is_read2][1]:
                        dummy,sc,xx,c1,c2,a,b,k1,k2,st,is_reverse = row
                        print("\t".join(map(str,[dummy,sc,xx,c1,c2,a,b,k1,k2] + [k1 in seq,k2 in seq,"{}:{}".format(seq,qal),st,is_reverse])))
                           # second_pass_reads[template_id,aln.is_secondary][1].append( ["gP:",sc,int(xx),c1,c2,c1[-1],c2[-1],k1,k2] )

    if args.debug: print("nreads:",nh)

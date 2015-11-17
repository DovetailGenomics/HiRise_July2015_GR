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
import random
#import chicago_edge_scores as ces
import pysam
#from random import random
from bisect import bisect
from bamtags import BamTags


cum_weights=[]
total_weight=0.0

debug=True

def weighted_choice(l,weighting):
    global cum_weights
    global total_weight
    if len(cum_weights)==0:
        total = 0
        for i in l:
            total += weighting[i]
            cum_weights.append(total)
        total_weight=total
    x = random.random() * total_weight
    i = bisect(cum_weights, x)
    if debug: print(i,l[i])
    return l[i]

def get_nlinks_dict(c1,bams,name2index,args,ll={}):

    count={}
#    index2=name2index[c2]
    for sam in bams:
        print("{}:{}-{}".format(c1,0,int(ll[c1])-1))
        for aln in sam.fetch(reference=c1):
            if (not aln.is_duplicate) and (aln.mapq >= args.mapq) and (BamTags.mate_mapq(aln) >= args.mapq):
                count[aln.rnext]=count.get(aln.rnext,0)+1
    return count

def get_nlinks(c1,c2,ll,bams,name2index,args):
    l1=ll[c1]
    l2=ll[c2]
    if l2<l1:
        x = c1
        c1= c2
        c2= x

    count=0
    index2=name2index[c2]
    for sam in bams:
        for aln in sam.fetch(region=c1):
            if (aln.rnext==index2) and (not aln.is_duplicate) and (aln.mapq >= args.mapq) and (BamTags.mate_mapq(aln) >= args.mapq):
                count+=1

    return count

def get_nlinks_spanning(c1,xx,ll,bams,name2index,args,gap):
    l1=ll[c1]

    count=0
#    index2=name2index[c2]
    for sam in bams:
        for aln in sam.fetch(region=c1):
            if (aln.rnext==aln.tid) and (not aln.is_duplicate) and (aln.mapq >= args.mapq) and (BamTags.mate_mapq(aln) >= args.mapq) and (aln.pos<xx-old_div(gap,2)) and (aln.pnext>xx+old_div(gap,2)):
                count+=1
    
    return count

def main():

    import argparse
    parser = argparse.ArgumentParser()
#    parser.add_argument('-l','--links')
    parser.add_argument('-P','--param')
    parser.add_argument('-o','--outfile')
    parser.add_argument('-b','--bamfiles',required=True)
    parser.add_argument('-N','--ncontigs',default=25,type=int)
    parser.add_argument('-q','--mapq',default=10,type=int)
    parser.add_argument('-G','--genomesize',default=3e9,type=float)
#    parser.add_argument('-n','--ncontigs',type=int,required=True)
    parser.add_argument('-d','--debug',default=False ,  action='store_true')
    parser.add_argument('-p','--progress',default=False ,  action='store_true')
    parser.add_argument('--seed',required=False,type=int,default=1, help="Seed for random number generation, use -1 for no seed")
    args = parser.parse_args()
    #ces.debug=args.debug
    if args.seed != -1 :
      random.seed(args.seed)

    if args.debug:
        args.progress=True

    
    oname={}
    bams=[]
    for bamfile in args.bamfiles.split(','):
        bams.append( pysam.Samfile(bamfile,"rb") )
        oname[bams[-1]]=bamfile
    h=bams[0].header
    seqs=h['SQ']

    slen=[ s['LN'] for s in seqs ]
    snam=[ s['SN'] for s in seqs ]
    contigs = snam
    ncontigs = len(snam)

    ll={}
    name2index={}
    G=0.0
    for i in range(len(snam)):
        if args.debug: print("n:",len(snam),i,snam[i],slen[i])
        ll[snam[i]]=float(slen[i])
        G+=ll[snam[i]]
        name2index[snam[i]]=i
    if args.debug: print("done")
#    G=args.genomesize

    def area(l1,l2):
        return (old_div(l1,1000.0))*(old_div(l2,1000.0))

    rate=0.0
    totaln=0.0
    n_pairs=0
    gap=300.0

    while n_pairs < args.ncontigs:
        n_pairs+=1
        c1 = weighted_choice(contigs,ll)
        if args.debug: print(c1)
#        c2 = weighted_choice(contigs,ll)
        totalnd_dict=get_nlinks_dict(c1,bams,name2index,args,ll)   #0.0+max(links.get((c1,c2),0), links.get((c2,c1),0))
        rate_list=[]

        for i in range(ncontigs):
            if not i == name2index[c1]:
                rated = area(ll[c1],slen[i])
                totalnd = totalnd_dict.get(i,0)
                rate_list.append( ( totalnd , rated ) )

        rate_list.sort()
        for i in range(len(rate_list)-5):
            totaln  += rate_list[i][0]
            rate    += rate_list[i][1]
        print("#inter:",1.0e-6*totaln/rate,G*G*1.0e-6*totaln/rate,G*1.0e-6*totaln/rate,G)
            
    print("final:",1.0e-6*totaln/rate,G*G*1.0e-6*totaln/rate,G*1.0e-6*totaln/rate,G)

    if args.outfile:
        f=open(args.outfile,"wt")
        f.write("{}\n".format(str({ 'G':G, 'Nn':G*G*1.0e-6*totaln/rate })))
        f.close()

if __name__=="__main__":
    main()

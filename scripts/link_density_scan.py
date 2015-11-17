#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
from hirise_assembly import HiriseAssembly
import struct
import hashlib
    
def process_links(links,w=1000,min_links=2,max_others=5,segments=False,debug=False):
    filtered_links={}
    blacklist={}
    l=[]
    for c1,c2 in list(sorted(links.keys())):
        #print c1,c2,links[c1,c2]
        filtered_links[c1,c2]=[]
        blacklist[c1,c2]=[]
        for x1,x2 in links[c1,c2]:
            l.append( (x1,c2,x2) )
    l.sort()
    #print l
    i=0
    j=0
    buff_stats={}
    bad_ones={}
    bad_steps=[]
    while i < len(l):
        #print j,i,l[j],l[i],ll[l[j][1]],ll[c1] #,filtered_links[c1,l[j][1]]
        buff_stats[l[i][1]] = buff_stats.get(l[i][1],0)+1
        while (l[i][0]-l[j][0])>w :
            buff_stats[l[j][1]] -= 1
            if buff_stats[l[j][1]]==0:
                del buff_stats[l[j][1]]
#            filtered_links[c1,l[j][1]].append( (l[j][0],l[j][2]) )
            j+=1
        ncs = sum( [ 1 for k in list(sorted(buff_stats.keys())) if buff_stats[k]>=min_links ] )
        i+=1
        if ncs>max_others:
#            for k in range(j,i):
#                bad_ones[k]=1
            bad_steps.append((j,+1))
            bad_steps.append((min(i,len(l)-1),-1))
            if debug: print("#bad:", j,i)

    bad_steps.sort()

    x=0
    last=0
    rs=0
    state="out"
    startx=0


    for i in range(len(bad_steps)):
        #print(i,bad_steps[i],len(l))
        y =bad_steps[i][0]
        x=l[y][0]
        if rs>0 and state=="out":
            state="in"
            startx=x
        rs+=bad_steps[i][1]        
        if rs<=0 and state=="in":
            if segments: 
#                segments.write("{} {} {}\n".format(c1,startx,x))
                print(c1,startx,x,"promiscuous",file=segments,sep="\t")
            state="out"



if __name__=="__main__":
     import sys
     import argparse

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug',default=False,action="store_true",help="Turn on debugging ouput")
     parser.add_argument('-C','--nchunks' ,default=1,type=int,help="Number of chunks.")
     parser.add_argument('-w','--window' ,   default=1000,type=int,help="Window size.")
     parser.add_argument('--mask'   ,default=[], action="append",help="Name of a file containing regions of the input assembly to mask out")
     parser.add_argument('-m','--min_links' ,default=2,type=int,help="Min links to another contig to count.")
     parser.add_argument('-M','--max_others',default=4,type=int,help="Max other contigs in the window.")
     parser.add_argument('-c','--chunk',default=0,  type=int,help="This chunk.")
     parser.add_argument('-q','--mapq',default=50,  type=int,help="Min map quality score.")
     parser.add_argument('-i','--infile',default=False,help="Filename for serialised assembly input file.")
     parser.add_argument('-o','--outfile',default=False,help="Filename for writing a list of segments on the raw contigs to mask for being promiscuous in linking.")
# -m 2 -w 1000 -M $( cat {input.threshold} ) 
     args = parser.parse_args()

     if args.infile:
          asf = HiriseAssembly()
          asf.load_assembly(args.infile)
     
     for segments_file in args.mask:
        asf.add_mask_regions(filename=segments_file)
        asf.merge_masked_regions()

     if args.outfile:
          of=open(args.outfile,"wt")
     else:
          of=sys.stdout
     for ocontig in asf.ocontigs_iter():
          chunk = struct.unpack("<L", hashlib.md5(ocontig.encode("utf-8")).digest()[:4])[0]%args.nchunks
          if not chunk == args.chunk: continue
          links={}
          asf.get_links([ocontig],skipI=True,mapq=args.mapq,links=links,contigs=False,raw=True)
#          print(ocontig,chunk,links)
          process_links(links,args.window,args.min_links,args.max_others,segments=of)
     of.close()

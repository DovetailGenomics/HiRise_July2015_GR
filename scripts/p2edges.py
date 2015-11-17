#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
from builtins import map
import sys
import string
import sys
import argparse
import idGen as idGen


if __name__=="__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d','--debug',     default=False, action="store_true", help="Turn on debugging ouput")
    parser.add_argument('-p','--progress',  default=False, action="store_true", help="Print progress info")
    parser.add_argument('-b','--blacklist', default=False,                      help="File containing edges to exclude.")
    parser.add_argument('-B','--breaks',    default=False,                      help="File containing breaks to add, in global coordinates.")
    
    args = parser.parse_args()
    if args.progress: print("#",args)

    breaks={}
    if args.breaks:
        fh = open(args.breaks)
        while True:
            l=fh.readline()
            if not l: break
            c=l.strip().split()
            breaks[c[0]]=breaks.get(c[0],[])+ [ c ]

    for b in breaks.keys():
        print(b)

    bl = {}
    if args.blacklist:
        fh = open(args.blacklist)
        while True:
            l=fh.readline()
            if not l: break
            c=l.strip().split()
            bl[c[3],c[4]]=1
            bl[c[4],c[3]]=1

    scaffold_hash={}
    last_c=-1
    last_x=0
    last_s=0
    cl=0
    contign=0
    scaffold_rename={}
    while True:
        l=sys.stdin.readline()
        if not l: break
        if not l[:2]=="p:": continue
        c=l.strip().split()
        scaffold_raw,end,x = c[1],c[2],float(c[3])
        scaffold = scaffold_rename.get(scaffold_raw,scaffold_raw)

        if not scaffold == last_s:
            #outfasta.next(name_prefix+"_"+str(scaffold))
            last_x=0
            contign=0
        contig,p=end[:-2],end[-1:]

        if not contig in scaffold_hash: 
            print("scaffold:",contig,scaffold)
            scaffold_hash[contig]=scaffold
        if (not contig == last_c) and last_s==scaffold:
            if not (contig,last_c) in bl:
                print("\t".join(map(str,["#edge:",last_end,end,{'length':int(x-last_x),'contig':False}])))

        elif last_c==contig:

            if args.breaks and is_broken(contig,breaks):

                new_scaffold = scaffold_raw + idGen.id()
                scaffold_rename[scaffold_raw]=new_scaffold
                if p=="3" : ## contig on the + strand

                    
#                      
#                    new_3prime_end = end
                    
                    
                    print("\t".join(map(str,["#edge:",last_end,end,{'length':int(x-last_x),'contig':True}])))
                    print("\t".join(map(str,["#edge:",last_end,end,{'length':int(x-last_x),'contig':True}])))

                    end     = new_contig + ".3"
                    contig  = new_contig
                    scaffold= new_scaffold
                else:       ## contig on the - strand
                    print("\t".join(map(str,["#edge:",last_end,end,{'length':int(x-last_x),'contig':True}])))
                    print("\t".join(map(str,["#edge:",last_end,end,{'length':int(x-last_x),'contig':True}])))                    
            else:
                print("\t".join(map(str,["#edge:",last_end,end,{'length':int(x-last_x),'contig':True}])))

        last_c=contig
        contign+=1
        last_x=x
        last_s=scaffold
        last_end=end


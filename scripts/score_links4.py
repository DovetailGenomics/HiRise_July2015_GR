#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
from builtins import str
from builtins import map
import sys
import networkx as nx
import chicago_edge_scores as ces
import numpy as np
from scipy.stats import poisson
import math

#G=3.0e9
N=100000000.0
pn=0.3
G=3000000000.0

import math
if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()
#    parser.add_argument('-t','--threshold',default=0.0 ,  type=float)
#    parser.add_argument('-b','--besthits')
    parser.add_argument('-d','--debug',default=False ,  action='store_true')
#    parser.add_argument('-P','--plot',default=False ,  action='store_true')
    parser.add_argument('-M','--set_insert_size_dist_fit_params',default=False )
#    parser.add_argument('-G','--maxTrueGap',type=int,default=False)
    parser.add_argument('-S','--scoreDelta',type=float,default=2.0)
    parser.add_argument('-p','--pvalue',type=float,default=0.000001)
    parser.add_argument('-E','--endwindow',type=float,default=False,help="Ignore links where either read is burried more than endwindow bases from the end of its mapped contig.")
    parser.add_argument('-N','--maxN',type=int,default=False)
    parser.add_argument('-L','--minlen',type=int,default=1000)

    args = parser.parse_args()
    ces.debug=args.debug
    fmodel=open( args.set_insert_size_dist_fit_params )
    contents = fmodel.read()
    try:
        fit_params=eval(contents)
    except:
        "couldn't deal with option", args.param
    fmodel.close
    ces.set_exp_insert_size_dist_fit_params(fit_params)


#    besthit={}
#    if args.besthits:
#        besthit={}
#        if args.besthits:
#            f = open(args.besthits)
#            while True:
#                l = f.readline()
#                if not l: break

#                if not l[:5]=="best:": continue
#                c=l.strip().split()
#                besthit[c[1]]=c[2:]
    #            print c[1],besthit[c[1]]
#            f.close()
#    if args.progress: print("#Done reading besthits")

    def ppf_cached(y,cache={}):
        x=round(y,4)
        if x==0:
            return(poisson.ppf(0.99999999,y))
        if x in cache: return cache[x]
        ppf=poisson.ppf(0.99999999,x)
        if np.isnan(ppf) or ppf==np.inf:
            print("wtf:",y,x)
        cache[x]=max(ppf,1)
        return ppf

    def ppf_mtc(y,N,cache={}):
        pp=(1.0-(args.pvalue/N))
        if pp==1.0:
            pp=1.0-np.finfo(float).resolution
        x=round(y,4)
        if x==0:
            ppf=poisson.ppf(pp,y)
            if np.isnan(ppf) or np.isinf(ppf):
                print("wtf:",y,x,ppf,(1.0-(args.pvalue/N)),N)
            return ppf
        if (x) in cache: return cache[x]
        ppf=poisson.ppf(pp,x)
        cache[x]=max(ppf,1)
        return ppf

    n_done=0
    G2=(ces.model.G*ces.model.G)
    while (not args.maxN) or n_done<args.maxN:
        l=sys.stdin.readline()
        if not l: break
        if l[0]=="#": continue
        c = l.strip().split("\t")
        s1,s2,l1,l2,n = c[0],c[1],int(c[2]),int(c[3]),int(c[4])
        if l1<args.minlen : continue
        if l2<args.minlen : continue
        if args.endwindow and ((l1>args.endwindow*2) or (l2>args.endwindow*2)) :
            links = eval(  " ".join(c[5:]) )
            n=0
            for x,y in links:
                if not (x<args.endwindow or (l1-x)<args.endwindow): continue
                if not (y<args.endwindow or (l2-y)<args.endwindow): continue
                n+=1
            l1 = min(l1,2*args.endwindow)
            l2 = min(l2,2*args.endwindow)
        n_done+=1
        n_bar0= ces.model.N*ces.model.pn*l1*l2*2/G2
#        n=len(links)
#        ppf=ppf_cached(n_bar0)
        N=int((G2/(l1*l2)))
        ppf=ppf_mtc(n_bar0,N)
        if (n-ppf)>=args.scoreDelta:
            print("\t".join(map(str,[s1,s2,n,ppf,round(n_bar0,4),l1,l2,N])))
 



#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
from builtins import range
import networkx as nx
import sys

import greedy_chicagoan2

def edges2oo(f):
    g=nx.Graph()
    ll={}
    scaffolds={}
    scaffold={}
    while True:
        l = f.readline()
        if not l: break
        if l[:6]=="#edge:": #continue
            c = l.strip().split()
            v=eval(" ".join(c[3:]))
            if v['contig']: ll[c[1][:-2]]=v['length']
            g.add_edge(c[1],c[2],v)
        if l[:3]=="cc:": #continue
            c = l.strip().split()
            scn=int(c[1])
            scl=eval(" ".join(c[3:]))
            scaffolds[scn]=scl
            for s in scl:
               scaffold[s]=scn 

    contigs=[]
    strand={}
    scaffold={}
    coords={}
    facing={}
    ccn=1
    for c in nx.connected_components(g):
#        print "#",len(c)
        ends=[]
        for cc in c:
            scaffold[cc]=ccn
            if g.degree(cc)==1:
                ends.append(cc)
        ccn+=1
        order = list(nx.dfs_preorder_nodes(g,ends[0]))

#        def traverse_and_layout(n,coords,facing,x,s,og,max_disp=False):
#    """Traverse the nodes in og, from node n.  Give node n position x.  s==1 for increasing coords, -1 for decreasing.  store in coords and facing the position and 'side' of each end respectively.  
#    Stop traversing if you hit max_disp (optional)."""
        greedy_chicagoan2.traverse_and_layout(ends[0],coords,facing,0,1,g)
        order1=[]
        for i in range(0,len(order),2):
            print(order)
            if not order[i][:-2]==order[i+1][:-2]:
                print("wtf", i,order)
                exit(0)
            if order[i][-1]=="5": 
                strand[order[i][:-2]]="+"
            else:
                strand[order[i][:-2]]="-"
            order1.append(order[i][:-2])
        contigs.append(order1)
    return(contigs,strand,ll,scaffolds,scaffold,coords)

if __name__=="__main__":

    import argparse
    parser = argparse.ArgumentParser()

    #parser.add_argument('-i','--input',default="gj4.out")
    parser.add_argument('-l','--links')
    parser.add_argument('-w','--window',type=int,default=3)

    args = parser.parse_args()
    print("#",args)

    contigs,strand,ll,scaffolds,scaffold,coords = edges2oo(sys.stdin)
    for cl in contigs:
        sl=max( [ coords[c+e] for c in cl for e in (".5",".3")] )
        cl.sort(key=lambda x: coords[x+".5"])
        for c in cl:
            if coords[c+".5"] < coords[c+".3"]:
                ends = (".5",".3")
            else:
                ends = (".3",".5")
                
            for e in ends:
                print("p:",scaffold[c+e],c+e,coords[c+e],"-",-1,sl,ll[c],1==2)
                #p: 1 Scaffold76818_1.3 790422 - -1 5832906 32362 False

#        print [strand[cc] for cc in c]

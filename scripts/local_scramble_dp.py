#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
from builtins import range
from builtins import object
import math
import itertools
import networkx as nx

class LSbitlinks(object):
    def __init__(self,w):
        n = 2**w
        self.n=n
        #print n
        bitlinks={}
        for d in range(1,w+1):
            for f in range(n):
                bitlinks[f,d]=[]
                for g in range(n):
                    #print d,f,g,bin(f),bin(g),bin(f>>d),bin((2**(w-d)-1)),bin(g&(2**(w-d)-1)),f>>d == g&(2**(w-d)-1)
                    if (f>>d) == g&(2**(w-d)-1):
                        bitlinks[f,d].append(g)
                #print d,f,bitlinks[f,d]
        self.bitlinks = bitlinks
        
class LSbacklinks(object):
    def __init__(self,w):

        g=nx.DiGraph()
        g0=nx.Graph()

        n=math.factorial(w)
        self.n = n
        perms = []
        self.perms=perms
        for p in itertools.permutations(list(range(w))):
            perms.append(p)

        for x in range(0,w)  :
            for y in range(x+1,w+1):
                d = y-x
                for i in range(n):
                    for j in range(n):
                        match=True
                        #print d,x,y,i,j,w-d,"####"
                        for k in range(w-d):
                            #print d,x,y,i,j,k,k+d
                            if not perms[i][k]+d ==  perms[j][k+d]:
                                match=False
                        if match:
                            g.add_edge((x,i),(y,j),weight=d*d )
                            g0.add_edge((x,i),(y,j),weight=d*d )
                            #print (x,i),(y,j),d,"z"

                            #        for e in g.edges():
                            #print "e:",e

        n_backlinks=0

        backlinks={}

        for i in range(n):
            backlinks[i]=[]
            nn = nx.dfs_preorder_nodes(g,(0,i))
            #print "#",str((0,i))#,list(nn)
            sg = g0.subgraph(nn)
            #print "#",list(sg.nodes())
            t = nx.minimum_spanning_tree(sg)
            for a,b in t.edges():
                if b<a:
                    c=b
                    b=a
                    a=c
                if a==(0,i):
                    #print "x:",b[0]-a[0],a,b,perms[b[1]],perms[a[1]]
                    n_backlinks+=1
                    backlinks[i].append( (b[0]-a[0],b[1]) )

        self.backlinks = backlinks
        self.n_backlinks=n_backlinks
        
if __name__=='__main__':
    import sys

    w=int(sys.argv[1])

    bl = LSbacklinks(w)
    
    #print "#",bl.n_backlinks,bl.n,float(bl.n_backlinks)/bl.n

    x="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"

    def sperm(p,s):
        r=""
        for pp in p:
            if pp<len(s):
                r+=s[pp]
        return r

    print(x)
    if False:
        for i in range(w,len(x)-w+1):
            for p in bl.perms:
                s1p = sperm(p,x[i:i+w])
                for o,p2 in bl.backlinks[p]:
                    s2 = x[i-o:i]
                    print("."*(i-o) + sperm(p2,s2) + s1p +"."*(len(x)-(i-o)-w-o),"%")
    else:
        for i in range(w,len(x)-w+1):
            for pi in range(bl.n):
                s1p = sperm(bl.perms[pi],x[i:i+w])
                for o,p2 in bl.backlinks[pi]:
                    s2 = x[i-o:i]
                    print("."*(i-o) + sperm(bl.perms[p2],s2) + s1p +"."*(len(x)-(i-o)-w-o),"%")



    fl = LSbitlinks(w)






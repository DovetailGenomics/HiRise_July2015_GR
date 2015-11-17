#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import map
from builtins import object
import sys
import networkx as nx
import greedy_chicagoan as gs
import math
import random
default_gapsize=100

def same_component(s1,s2,ccd):
#    print "same component?",ccd[s1],ccd[s2]
    cs1=ccd[s1]
    while cs1 in ccd: cs1=ccd[cs1]
    cs2=ccd[s2]
    while cs2 in ccd: cs2=ccd[cs2]
    if cs1==cs2: return True
    return False

def merge_components(s1,s2,ccd):
    cs1=ccd[s1]
    while cs1 in ccd: cs1=ccd[cs1]
    cs2=ccd[s2]
    while cs2 in ccd: cs2=ccd[cs2]

    if not cs1==cs2:
        ccd[cs1]=cs2
    
#        print cs2,cs1,"ccd[{}]={}".format(cs1,cs2)


class ScaffoldEdit(object):
    def __init__(self,repr):
        if type(repr)==tuple:
            self.score=repr[0]
            self.breaks=repr[1]
            self.joins=repr[2]
        elif type(repr)==dict:
            self.score=repr['score']
            self.breaks=repr.get('breaks',[])
            self.joins=repr['joins']

    def is_valid(self,links,ccd):
#        print "#valid?",self
        freed=[]
        for a,b in self.breaks:
            if not (a,b) in links: return  False            # False
            freed.append(a)
            freed.append(b)

        if len(self.joins)==0: return True
#        print "#"
#        print self.joins
        for a,b in self.joins:
#            print a,b,a in links, b in links, a in freed, b in freed
#            print b in freed
            if a in links and not a in freed: return False  # False
            if b in links and not b in freed: return False  #False
            if same_component(a,b,ccd): return False        #False
        return True

    def implement(self,links,ccd,g):
#        return
#        print "implement:",self
        for a,b in self.breaks:
            if (a,b) in links: del links[a,b]
            if (b,a) in links: del links[b,a]
            del links[a]
            del links[b]
            if g.has_edge(a,b):
                g.remove_edge(a,b)
            else:
                print("how come there's not edge ?",a,b)
                raise Exception 
#            if (a,b) in g: g.remove_edge(a,b)

        for a,b in self.joins:
            links[a,b]=1
            links[b,a]=1
            links[a]=b
            links[b]=a
#            g.add_edge(a,b,weight=self.score)
            g.add_edge(a,b, {'length': default_gapsize, 'contig': False} )
            merge_components(a,b,ccd)
            if not (a,b) in links: isv=False

    def __repr__(self):
        return "\t".join(map(str,[ self.score, self.joins, self.breaks ]))

def update_end_distance(end_distance,n,g):
    x=0
    q=[n]
    seen={}
    last=False
    while len(q)>0:
        m=q.pop(0)
        if last:
            try:
                x+=g[last][m]['length']
            except Exception as e:
                print("wtf?",last,m,x,n)
                print(math.log(-1.0))
        last=m
        if (not m in end_distance) or end_distance[m]>x: end_distance[m]=x
        
        seen[m]=True
        if len(g.neighbors(m))>2:
            print("topology error:",m,list(g.neighbors(m)))
            print(math.log(-1.0))
        for l in g.neighbors(m):
            if not l in seen:
                q.append(l)

    return x

#    print end_distance

#def llr(e1,e2):
#    return 1.0

L=200000.0

if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-d','--debug',default=False,action='store_true')
    parser.add_argument('-p','--progress',default=False,action='store_true')
#    parser.add_argument('-L','--links')
    parser.add_argument('-s','--scaffolds')
#    parser.add_argument('-S','--alreadyDone')
    parser.add_argument('-b','--besthits')
    parser.add_argument('-l','--lengths')
    parser.add_argument('-E','--edgefile')
    parser.add_argument('-F','--filter')
    parser.add_argument('-N','--name')
    parser.add_argument('-m','--minscore' ,default=5.0,type=float)
    parser.add_argument('--seed',required=False,type=int,default=1, help="Seed for random number generation, use -1 for no seed")
#    parser.add_argument('-K','--slices'   ,default=1,type=int)
#    parser.add_argument('-k','--slice'    ,default=0,type=int)

    args = parser.parse_args()
    if args.seed != -1 :
      random.seed(args.seed)
    if args.debug:
        args.progress=True

    if args.progress: log( str(args) )

    name_prefix=""
    if args.name:
        name_prefix=args.name
    else:
        import idGen as idGen
        name_prefix="Scf" + idGen.id()

    ll={}
    if args.lengths:
        f = open(args.lengths)
        while True:
            l = f.readline()
            if not l: break
            if l[0]=="#": continue

            c=l.strip().split()
            ll[c[0]]=int(c[1])
        f.close()

    besthit={}
    if args.besthits:
#        besthit={}
        if args.besthits:
            f = open(args.besthits)
            while True:
                l = f.readline()
                if not l: break

                if not l[:5]=="best:": continue
                c=l.strip().split()
                besthit[c[1]]=c[2:]
    #            print c[1],besthit[c[1]]
            f.close()
    if args.progress: print("#Done reading besthits")

    linked={}
    g=nx.Graph()
    if args.scaffolds:
        f=open(args.scaffolds)
        while True:
            l=f.readline()
            if not l: break
            c=l.strip().split()
            if c[0]=="#edge:":
                at=eval(" ".join(c[3:]))
                a,b=c[1],c[2]
                g.add_edge(a,b,at)
                if not at['contig']:
                    linked[a]=1
                    linked[b]=1
                    linked[a,b]=1
                    linked[b,a]=1
#                print "#add edge",c[1],c[2],eval(" ".join(c[3:]))
    sys.stdout.flush()
    sc=1
    scaffold={}
    for c in nx.connected_components(g):
        for cc in c:
            scaffold[cc]=sc
            scaffold[cc[:-2]]=sc
        sc+=1


#    scaffold_pairs_tested={}
#Scaffold50016_1 Scaffold40593_1 ['Scaffold77744_1.5', 'Scaffold246520_1.5'] ['Scaffold111955_1.3', 'Scaffold216064_1.3'] 1141 455 1 15 1

#    joins_g=nx.Graph()

    moves=[]

    #={}
    while True:
        l=sys.stdin.readline()
        if not l: break
#        print l
#        print "\""+l[:10]+"\""
        if l[:10]=="link score": 
            c=l.strip().split()
            x=max(list(map(float,c[4:])))
            if x > args.minscore:
                moves.append( ScaffoldEdit({ 'score':x , 'joins': ((c[2],c[3]),)  }) )

        elif l[:7]=="interc:": 
            c=l.strip().split()
            x=eval(" ".join(c[1:]))
            if x[0] > args.minscore:
                moves.append( ScaffoldEdit({ 'score':x[0] , 'breaks': x[2]  , 'joins': x[1]  }) )

    moves.sort(key=lambda x: x.score,reverse=True)

    cnx=1
    ccd={}
    for c in nx.connected_components(g):
        for cc in c:
            ccd[cc]=cnx
        cnx+=1

    for m in moves:
        print(m,m.is_valid(linked,ccd))
        if m.is_valid(linked,ccd)==True:
            m.implement(linked,ccd,g)

#    exit(0)

    end_distance={}
    sn=1
    for sg in nx.connected_component_subgraphs(g):
        ends=[]
        bh_stats={}
        for n in sg.nodes():
            if sg.degree(n)==1:
                ends.append(n)

        if len(ends)==0: 
            print("why no ends?", sn) 
            sn+=1
            continue
        maxx=update_end_distance(end_distance,ends[0],sg)

# ['34329.0', '3', '+', '71834554', '71853152', '1', '1', '18598']

        t=0
        gap_len=0
        for s1,s2 in sg.edges():
            t+=sg[s1][s2]['length']
            if not sg[s1][s2]['contig']: gap_len += sg[s1][s2]['length']
            print("#",sn,s1,s2,sg[s1][s2]['length'],t)
        print(t,n,"slen",gap_len,t-gap_len)

        node_tlen=0
        nodes_list = list(sg.nodes())
        nodes_list.sort(key=lambda x: end_distance[x])
        for n in nodes_list:
            base_name,end_id=n[:-2],n[-1:]
            if end_id=="5": node_tlen+= ll[base_name]
            bh = besthit.get(base_name,False)
            x=-1
            chr="-"
            if bh:
                chr=bh[1]
                if bh[2]=="+":
                    if n[-1:]=="5":
                        x=int(bh[3])
                    else:
                        x=int(bh[4])
                if bh[2]=="-":
                    if n[-1:]=="5":
                        x=int(bh[4])
                    else:
                        x=int(bh[3])
                    
            print("p:",sn,n,end_distance[n],chr,x,t,ll[base_name],bh)
        print(node_tlen,"node_tlen")

        sn+=1

    exit(0)

    exit(0)


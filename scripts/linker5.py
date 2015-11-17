#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
from builtins import str
from builtins import map
from builtins import range
import sys
import networkx as nx
import greedy_chicagoan2 as gs
import glob
import random




def update_end_distance(end_distance,n,g):
    x=0
    q=[n]
    seen={}
    last=False
    while len(q)>0:
        m=q.pop(0)
        if last:
            x+=g[last][m]['length']
        last=m
        if (not m in end_distance) or end_distance[m]>x: end_distance[m]=x
        
        seen[m]=True
        for l in g.neighbors(m):
            if not l in seen:
                q.append(l)

#    print end_distance

#def llr(e1,e2):
#    return 1.0

L=200000.0

def test_interc(c1,c2,join_options,graph,ends,linked,slen,max_interc_len):
    os1=oscaffold.get(c1,-1)
    #os1_mine = (args.fromS <= os1) and (os1 < args.endS)
    #if not os1_mine: return([-100])
    #print "test interc:",c1,c2,oscaffold.get(c1,0),oscaffold.get(c2,0) ,ends[scaffold[c1]],ends[scaffold[c2]]#,c2_end in linked, c1_end in linked #,linked.get(c1),linked.get(c2),len(linked.keys())
    print("test interc:", c1,c2,scaffold[c1],scaffold[c2],ends[scaffold[c1]],ends[scaffold[c2]])
    if slen[scaffold[c1]] < max_interc_len: 
        for free_end in ends[scaffold[c1]]:
            for c2_end in [c2+".5",c2+".3"]:
#                print "c2_end in linked:",c2_end,c2_end in linked
                if c2_end in linked:
                    print(1,free_end,c2_end)
                    try:
                        gs.test_interc_option( c2_end, linked[c2_end], free_end, join_options, graph )
                    except Exception as e:
                        print("c2_end=",c2_end,"--",linked[c2_end],free_end)
                        print(e)
                        raise Exception('too connected 1')
    if slen[scaffold[c2]] < max_interc_len:
        for free_end in ends[scaffold[c2]]:
            for c1_end in [c1+".5",c1+".3"]:
#                print "c1_end in linked:",c1_end,c1_end in linked
                if c1_end in linked:
                    print(2,free_end,c1_end)
                    try:
                        gs.test_interc_option( c1_end, linked[c1_end], free_end, join_options, graph )
                    except Exception as e:
                        print("c1_end=",c1_end,"--",linked[c1_end],free_end)
                        print(e)
                        raise Exception('too connected 2')

    if len(join_options)>0:
        join_options.sort(reverse=True)
        print("join options:",join_options)
        return(join_options[0])
    return([-100.0])

if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-d','--debug',default=False,action='store_true')
    parser.add_argument('-p','--progress',default=False,action='store_true')
    parser.add_argument('-L','--links')
    parser.add_argument('-a','--chunk1',type=int)
    parser.add_argument('-b','--chunk2',type=int)
    parser.add_argument('-1','--first')
    parser.add_argument('-2','--second')
    parser.add_argument('-s','--scaffolds')
    parser.add_argument('-S','--alreadyDone')
    #parser.add_argument('-b','--besthits')
    parser.add_argument('-l','--max_interc_len',default=20000,type=int)
    parser.add_argument('-E','--edgefile')
    parser.add_argument('--test_intercs',default=False,action="store_true")
    parser.add_argument('-F','--filter')
    parser.add_argument('-K','--endS',default=100000,type=int)
    parser.add_argument('-k','--fromS' ,default=0,type=int)
    parser.add_argument('-N','--nchunks' ,default=32,type=int)
    parser.add_argument('-t','--threshold' ,default=0.10,type=float)
    parser.add_argument('-Z','--chunkfile')
    parser.add_argument('-M','--set_insert_size_dist_fit_params',default="3.85301461797326,1.42596694138494,1.38674994280385e-05,10940.8191219759,49855.7525034142,0.3,420110993")
    parser.add_argument('--seed',required=False,type=int,default=1, help="Seed for random number generation, use -1 for no seed")

    args = parser.parse_args()

    if args.seed != -1 :
        random.seed(args.seed)

    if args.debug:
        args.progress=True

    max_interc_len=args.max_interc_len

    if args.progress: print("#", str(args)) 

    print("#"+str(args))
    fit_params={}
    try:
        fit_params = eval(args.set_insert_size_dist_fit_params )
    except Exception as e:
        f=open( args.set_insert_size_dist_fit_params )
        contents = f.read()
        try:
            fit_params=eval(contents)
        except:
            "couldn't deal with option", args.param
        f.close
        gs.set_exp_insert_size_dist_fit_params(fit_params)

    ll={}

    oscaffold={}
    #oslen={}
    slen={}
    g=nx.Graph()
    linked={}

    scaffold_hashes={}
    my_contigs={}
    my_scaffolds={}
    scaffold={}
    import hashlib
    import struct

    while True:
        l=sys.stdin.readline()
        if not l: break
        c=l.strip().split()
        if c[0]=="#edge:":
            ht=eval(" ".join(c[3:]))
            g.add_edge(c[1],c[2],ht)
            if not ht['contig']:
                linked[c[1]]=c[2]
                linked[c[2]]=c[1]
            scaffold[c[1]]=scaffold[c[1][:-2]]
            scaffold[c[2]]=scaffold[c[2][:-2]]
            slen[scaffold[c[1]]]=slen.get(scaffold[c[1]],0)+ht['length']
            if ht['contig']: ll[c[1][:-2]]=ht['length']
        if c[0]=="scaffold:":
            contig,cscaffold=c[1],c[2]
            scaffold[contig]=cscaffold
            if not cscaffold in scaffold_hashes:
                h=struct.unpack("<L", hashlib.md5(cscaffold.encode("utf-8")).digest()[:4])[0]
                scaffold_hashes[cscaffold]=h%args.nchunks
#                print "scaffold_hash",c[1],scaffold_hashes[c[1]] 
                #print "scaffold_hash",cscaffold,scaffold_hashes[cscaffold],args.chunk1,args.chunk2
                if scaffold_hashes[cscaffold] in (args.chunk1, args.chunk2):
                    my_scaffolds[cscaffold]=1
            if scaffold_hashes[cscaffold] in (args.chunk1, args.chunk2):
                my_contigs[c[1]]=1
               # print "my contig",c[1],cscaffold
            

    linkfiles = list(glob.glob(args.links))
    random.shuffle(linkfiles)
    print(linkfiles)

    links={}
    links_interc={}
    nlinks=0
    scaffold_pair_links={}
    pairs_to_test={}
    inter_chunk_pairs={}
    inter_chunk_scaffolds={}
    for lfilename in linkfiles:
        print("#",lfilename)
        f = open(lfilename,"rt")
        while True:
            l = f.readline()
            if not l: break
            if l[0]=="#": continue
            c=l.strip().split()
            c1,c2=c[0],c[1]
            if not (c1 in my_contigs and c2 in my_contigs) : continue
            if not c1 in scaffold: continue
            if not c2 in scaffold: continue
            s1,s2=scaffold[c1],scaffold[c2]
            if s1==s2: continue
            if not((scaffold_hashes[s1]==args.chunk1 and scaffold_hashes[s2]==args.chunk2) or (scaffold_hashes[s1]==args.chunk2 and scaffold_hashes[s2]==args.chunk1)): continue
            links[c1,c2]=eval(" ".join(c[5:]))
            nlinks+=1
            if nlinks%1000==1: print("#nlinks=",nlinks,c1,c2,scaffold[c1],scaffold[c2])
            if not c1 in ll: ll[c1]=int(c[2])
            if not c2 in ll: ll[c2]=int(c[3])
            sp=tuple(sorted([scaffold[c1],scaffold[c2]]))
            scaffold_pair_links[sp]=scaffold_pair_links.get(sp,0)+1
            if not sp in inter_chunk_scaffolds:
                inter_chunk_pairs[c1,c2]=1
                inter_chunk_scaffolds[sp]=1
        f.close()
    end_distance={}
    ends={}

    for n in g.nodes():
        if g.degree(n)==1:# and n in oscaffold and oscaffold[n]<args.endS:
            update_end_distance(end_distance,n,g)
            ends[scaffold[n]] = ends.get(scaffold[n],[])+[n]
#x            print n,scaffold[n],ends[scaffold[n]]

    print("#done setting edge distances")
    sys.stdout.flush()

    max_interc_len = 20000
    gaps_list=(1000,)
    scaffold_pairs_tested={}
    link_scores={}
    gs.ll=ll
    gs.links=links
    pairs_to_test=list(inter_chunk_pairs.keys())
#    if args.first:
#        pairs_to_test = [(args.first,args.second)]
#    for c1,c2 in links.keys():
#    for c1,c2 in [("Scaffold68143_1","Scaffold42944_1"),("Scaffold232867_1","Scaffold519668_1"),("Scaffold82730_1","Scaffold59156_1")]:
    for c1,c2 in pairs_to_test:
        nlinks_used={}
        if c1 in scaffold and c2 in scaffold and scaffold[c1]==scaffold[c2]:   continue
        #if c1+".5" in nx.node_connected_component(g,c2+".5" ): print "wtf? 1 ",c1,c2

        #if oscaffold[c1]==oscaffold[c2]: print "wtf? 3 ",c1,c2,oscaffold[c1],oscaffold[c2]
#        if not scaffold[c1]%args.slices == args.slice : continue
#        scaffold_pairs_tested[c2,c1]=True
        #print oslen[oscaffold.get(c1,c1)]
        header="\t".join(map(str,[ c1,c2,ends[scaffold[c1]],ends[scaffold[c2]],scaffold[c1],scaffold[c2],oscaffold.get(c1),oscaffold.get(c2),slen[scaffold[c1]],slen[scaffold[c2]]])) #,scaffold[c1]%args.slices,scaffold[c2]%args.slices,args.slice
#        sys.stdout.flush()
        linkscore=[]
        hit=False

        if not ((scaffold[c1],scaffold[c2]) in scaffold_pairs_tested): 
            for e1 in ends[scaffold[c1]]:
                for e2 in ends[scaffold[c2]]:
    #                link_scores[e1,e2]=[ gs.link_test(g,e1,e2,gap) for gap in (1, 1000,5000,10000, 20000, 35000, 50000, 80000, 100000,150000,200000,500000 ) ]
                    link_scores[e1,e2]               =[ gs.link_test(g,e1,e2,gap) for gap in gaps_list ]
                    if args.debug: nlinks_used[e1,e2]=[ gs.nlinks_used(g,e1,e2,gap) for gap in gaps_list ] 
                    if max( link_scores[e1,e2])>args.threshold: hit=True

            print(ends[scaffold[c1]], slen[scaffold[c1]],slen[scaffold[c2]], max( link_scores[e1,e2]), max(nlinks_used.get((e1,e2),[-1])))

        scaffold_pairs_tested[scaffold[c1],scaffold[c2]]=True
        scaffold_pairs_tested[scaffold[c2],scaffold[c1]]=True

        best_interc=""

        if (slen[scaffold[c1]]<=max_interc_len) or (slen[scaffold[c2]]<=max_interc_len):
            if args.test_intercs:
                join_options=[]
                best_interc=[-1]
                #            try:
                best_interc=test_interc(c1,c2,join_options,g,ends,linked,slen,max_interc_len)
                print("#best_interc:",best_interc)
    #            except Exception as e:
    #                print e
    #                print header
                if best_interc[0]>0:
                    hit=True
                else:
                    pass
    #                best_interc=""
        #        def joined(c,graph):

#        if args.test_intercs:
        if hit:
            print(header,"i:",best_interc)
            if best_interc:
                print("interc:",best_interc)
            for e1 in ends[scaffold[c1]]:
                for e2 in ends[scaffold[c2]]:
                    if max(link_scores[e1,e2])>-100:
                        print("link score",e1,e2,"\t".join(map(str,link_scores[e1,e2])))
                    if args.debug:
                        for ii in range(len(gaps_list)):
                            print("#nlink_debug:",e1,e2,ii,nlinks_used[e1,e2][ii],link_scores[e1,e2][ii]) 
            sys.stdout.flush()

#                            sc[x,y] = link_test(og,x,y)
#link score Scaffold68143_1.5 Scaffold42944_1.3 -12.0363618331   -8.50975484023  3.09050796551   13.1466146475   23.6521102192
#link score Scaffold232867_1.5 Scaffold519668_1.3 14.2053772843  17.5334920011   28.5300425318   38.1409186507   48.2455365211
#link score Scaffold82730_1.3 Scaffold59156_1.3 19.032925025     22.7060494622   34.76613982     45.137277127    55.5327333446
#link score Scaffold139910_1.5 Scaffold88540_1.5 18.4718988438   22.0553877336   33.8485414561   44.045490051    54.4456783858
#link score Scaffold264145_1.3 Scaffold163676_1.5 58.5394818429  61.4418869438   70.6603831852   77.9258584168   83.0945142366
#link score Scaffold48407_1.3 Scaffold136888_1.5 43.5898317791   46.6738612148   56.8244456982   65.4217034062   73.3441397183
#link score Scaffold113693_1.5 Scaffold61032_1.5 23.6894794909   27.0904366266   38.3018179122   47.9797590877   57.7542115065
#link score Scaffold125405_1.5 Scaffold158227_1.3 24.2942327515  27.8092174693   39.3813872363   49.3913743569   59.5690031451

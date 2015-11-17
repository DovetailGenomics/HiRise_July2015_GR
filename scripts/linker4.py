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

def test_interc(c1,c2,join_options,graph,ends,linked):
    os1=oscaffold[c1]
    os1_mine = (args.fromS <= os1) and (os1 < args.endS)
    if not os1_mine: return([-100])
    #        print "test interc",c1,c2,oscaffold[c1],oscaffold[c2],ends[oscaffold[c1]],ends[oscaffold[c2]],linked.get(c1),linked.get(c2),len(linked.keys())
    if False:
        for free_end in ends[scaffold[c1]]:
            for c2_end in [c2+".5",c2+".3"]:
                if c2_end in linked:
                    #                    print 1,free_end,c2_end
                    try:
                        gs.test_interc_option( c2_end, linked[c2_end], free_end, join_options, graph )
                    except Exception as e:
                        print("c2_end=",c2_end,"--",linked[c2_end],free_end)
                        print(e)
                        raise Exception('too connected 1')
    for free_end in ends[scaffold[c2]]:
        for c1_end in [c1+".5",c1+".3"]:
            if c1_end in linked:
                #                    print 2,free_end,c1_end
                try:
                    gs.test_interc_option( c1_end, linked[c1_end], free_end, join_options, graph )
                except Exception as e:
                    print("c1_end=",c1_end,"--",linked[c1_end],free_end)
                    print(e)
                    raise Exception('too connected 2')

    if len(join_options)>0:
        join_options.sort(reverse=True)
        #            print join_options
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
    parser.add_argument('-l','--lengths')
    parser.add_argument('-E','--edgefile')
    parser.add_argument('--test_intercs',default=False)
    parser.add_argument('-F','--filter')
    parser.add_argument('-K','--endS',default=100000,type=int)
    parser.add_argument('-k','--fromS' ,default=0,type=int)
    parser.add_argument('-t','--threshold' ,default=0.10,type=float)
    parser.add_argument('-Z','--chunkfile')
    parser.add_argument('--set_insert_size_dist_fit_params',default="3.85301461797326,1.42596694138494,1.38674994280385e-05,10940.8191219759,49855.7525034142,0.3,420110993")

    args = parser.parse_args()
    if args.debug:
        args.progress=True

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

    components=[]
    component={}
    chunk={}
    if args.chunkfile:
        f = open(args.chunkfile)
        while True:
            l = f.readline()
            if not l: break
            c = l.strip().split()
            if not c[0]=="c:": continue
            if int(c[2])==args.chunk1 or int(c[2])==args.chunk2:
                ccn = int(c[1])
                contigs = eval(" ".join(c[4:]))
                components.append(contigs)
                for cc in contigs:
                    component[cc]=ccn
                    chunk[cc]=int(c[2])

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

    oscaffold={}
    oslen={}
    slen={}
    g=nx.Graph()
    linked={}
    if args.scaffolds:
        f=open(args.scaffolds)
        while True:
            l=f.readline()
            if not l: break
            c=l.strip().split()
            if c[0]=="cc:":
                #print l.strip()
#                cc: 10 313 ['Scaffold98750_1', 'Scaffold2044
                scn=int(c[1])
                scs=eval(" ".join(c[3:]))
                for s in scs:
                    oscaffold[s]=scn
                    g.add_node(s+".3")
                    g.add_node(s+".5")
                    g.add_edge(s+".3",s+".5",{"contig":True, "length":ll[s]})
            if c[0]=="#edge:":
                ht=eval(" ".join(c[3:]))
                g.add_edge(c[1],c[2],ht)
                sb,e=c[1][:-2],c[1][-1:]
                
                oslen[oscaffold.get(sb,sb)] = oslen.get(oscaffold.get(sb,sb),0)+ht['length']
                if not ht['contig']:
                    linked[c[1]]=c[2]
                    linked[c[2]]=c[1]
                #print "#add edge",c[1],c[2],eval(" ".join(c[3:])),linked[c[1]],linked[c[2]]
        f.close()

        for n in g.nodes():
            s = n[:-2]
            g.add_edge(s+".3",s+".5",{"contig":True, "length":ll[s]})
            if not s in oscaffold:
                oscaffold[s]=0

#    for n in ["Scaffold27527_1","Scaffold40367_1","Scaffold17670_1"]:
#        print n,oscaffold[n] #, nx.node_connected_component(g,n)

#    print len(linked.keys())
    sys.stdout.flush()
    sc=1
    scaffold={}
    for c in nx.connected_components(g):
        for cc in c:
            scaffold[cc]=sc
            scaffold[cc[:-2]]=sc
        tl=0
        for e1,e2 in nx.dfs_edges(g,c[0]):
            tl+=g[e1][e2]['length']
        slen[sc]=tl
        sc+=1


    end_distance={}
    ends={}

    for n in g.nodes():
        if g.degree(n)==1:# and n in oscaffold and oscaffold[n]<args.endS:
            update_end_distance(end_distance,n,g)
            ends[scaffold[n]] = ends.get(scaffold[n],[])+[n]
#x            print n,scaffold[n],ends[scaffold[n]]

    print("#done setting edge distances")
    sys.stdout.flush()

    scaffold_pairs_tested={}
#Scaffold50016_1 Scaffold40593_1 ['Scaffold77744_1.5', 'Scaffold246520_1.5'] ['Scaffold111955_1.3', 'Scaffold216064_1.3'] 1141 455 1 15 1
    if args.alreadyDone:
        f=open(args.alreadyDone)
        while True:
            l=f.readline()
            if not l: break
            if l[0]=="#": continue
            if "link score"==l[:10]: continue
            c=l.strip().split()
            print("# skip",int(c[-5]),int(c[-4]))
            s1=scaffold[c[0]]
            s2=scaffold[c[1]]
            scaffold_pairs_tested[s1,s2]=True
            scaffold_pairs_tested[s2,s1]=True
        f.close()

    max_interc_len = 20000

    inter_chunk_pairs={}

    links={}
    links_interc={}
    if args.links:
        if args.links=="-":
            f = sys.stdin
        else:
            f = open(args.links)
        nlinks=0
        while True:
            l = f.readline()
            if not l: break
            if l[0]=="#": continue
#    parser.add_argument('-K','--endS',default=100000,type=int)
#    parser.add_argument('-k','--fromS' ,default=0,type=int)

            c=l.strip().split()
            c1,c2=c[0],c[1]
            if ( chunk.get(c1)==args.chunk1 and chunk.get(c2)==args.chunk2 ) or ( chunk.get(c2)==args.chunk1 and chunk.get(c1)==args.chunk2 ) :
                inter_chunk_pairs[c1,c2]=1
            os1=oscaffold.get(c1,-1)
            os2=oscaffold.get(c2,-1)            
            os1_mine = (args.fromS <= os1) and (os1 < args.endS)
            os2_mine = (args.fromS <= os2) and (os2 < args.endS)
#            if (os1_mine or os2_mine)  
            if (os1_mine or os2_mine)  \
                  and os1>=0 and os2>=0                                                                     \
                  and (c1 in scaffold and c2 in scaffold)                                                   \
                  and (not ( scaffold[c1],scaffold[c2] ) in scaffold_pairs_tested):
                links[c1,c2]=eval(" ".join(c[5:]))
                nlinks+=1
                if nlinks%10000==1: print("#nlinks=",nlinks,c1,c2,os1,os2,scaffold[c1],scaffold[c2],args.fromS,args.endS)
                
            else:
                pass
                #print c1,c2,ll[c1],ll[c2],eval(" ".join(c[5:]))
            #if oslen[oscaffold.get(c1,c1)] < max_interc_len or oslen[oscaffold.get(c1,c1)] < max_interc_len:
            #    links_interc[c1,c2] = eval(" ".join(c[5:])) 
        if args.links=="-":
            pass
        else:
            f.close()

    gaps_list=(1000,)

    link_scores={}
    gs.ll=ll
    gs.links=links
    pairs_to_test=list(inter_chunk_pairs.keys())
    if args.first:
        pairs_to_test = [(args.first,args.second)]
#    for c1,c2 in links.keys():
#    for c1,c2 in [("Scaffold68143_1","Scaffold42944_1"),("Scaffold232867_1","Scaffold519668_1"),("Scaffold82730_1","Scaffold59156_1")]:
    for c1,c2 in pairs_to_test:
        nlinks_used={}
        if c1 in scaffold and c2 in scaffold and scaffold[c1]==scaffold[c2]:   continue
        if (scaffold[c1],scaffold[c2]) in scaffold_pairs_tested: continue
        #if c1+".5" in nx.node_connected_component(g,c2+".5" ): print "wtf? 1 ",c1,c2

        #if oscaffold[c1]==oscaffold[c2]: print "wtf? 3 ",c1,c2,oscaffold[c1],oscaffold[c2]
#        if not scaffold[c1]%args.slices == args.slice : continue
        scaffold_pairs_tested[scaffold[c1],scaffold[c2]]=True
        scaffold_pairs_tested[scaffold[c2],scaffold[c1]]=True
#        scaffold_pairs_tested[c2,c1]=True
        #print oslen[oscaffold.get(c1,c1)]
        header="\t".join(map(str,[ c1,c2,ends[scaffold[c1]],ends[scaffold[c2]],scaffold[c1],scaffold[c2],oscaffold.get(c1),oscaffold.get(c2),slen[scaffold[c1]],slen[scaffold[c2]]])) #,scaffold[c1]%args.slices,scaffold[c2]%args.slices,args.slice
#        sys.stdout.flush()
        linkscore=[]
        hit=False
        for e1 in ends[scaffold[c1]]:
            for e2 in ends[scaffold[c2]]:
#                link_scores[e1,e2]=[ gs.link_test(g,e1,e2,gap) for gap in (1, 1000,5000,10000, 20000, 35000, 50000, 80000, 100000,150000,200000,500000 ) ]
                link_scores[e1,e2]               =[ gs.link_test(g,e1,e2,gap) for gap in gaps_list ]
                if args.debug: nlinks_used[e1,e2]=[ gs.nlinks_used(g,e1,e2,gap) for gap in gaps_list ] 
                if max( link_scores[e1,e2])>args.threshold: hit=True

        print(ends[scaffold[c1]], slen[scaffold[c1]],slen[scaffold[c2]], max( link_scores[e1,e2]), max(nlinks_used.get((e1,e2),[-1])))

        best_interc=""
        if args.test_intercs:
            join_options=[]
            best_interc=[-1]
            try:
                best_interc=test_interc(c1,c2,join_options,g,ends,linked)
            except Exception as e:
                print(e)
                print(header)
            if best_interc[0]>0:
                hit=True
            else:
                best_interc=""
    #        def joined(c,graph):

        if hit:
            print(header,"i:",best_interc)
            for e1 in ends[scaffold[c1]]:
                for e2 in ends[scaffold[c2]]:
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

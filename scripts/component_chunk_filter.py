#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import map
from builtins import str
from builtins import range
from past.utils import old_div
import sys
import networkx as nx
import chicago_edge_scores as ces
import random

def is_shaved_tail(G,shave_round,shaved_degree,shave_limit):
    leaves=[]
    r={}
    for n in sorted(G.nodes()):
        if shave_round.get(n,0)>shave_limit and shaved_degree.get(n,0)==1:
            leaves.append(n)
    for l in sorted(leaves):
        q=[l]
        while len(q)>0:
            n=q.pop()
            r[n]=True
            for nn in G.neighbors(n):
                if shave_round.get(nn,0)>shave_limit and shaved_degree.get(nn,0)<=2 and (not nn in q) and (not nn in r ):
                    q.append(nn)
    return r


def distance_to_nearest_branch(G,shave_round,shave_limit,trim_degree):
    branchpoints=[]
    r={}

    for n in G.nodes():
        #print "#",n, G.degree(n),G.degree(n)==1
        if shave_round.get(n,0)>shave_limit and trim_degree.get(n,0)>2:
            branchpoints.append(n)
#    print leaves

    round=1
#    print "#",leaves
    boundary=list(branchpoints)
    next=[]
    done=[]
    #r={}
    for b in boundary: r[b]=round
    while len(boundary)>0:
        round +=1
        next=[]
        for n in boundary:
            for nn in G.neighbors(n):
                if (not nn in boundary+done+next) :
                    next.append(nn)
#                    r[nn]=round
        for b in next: 
            r[b] = round
#            if len( [nn for nn in G.neighbors(b) if not r.has_key(nn)] )==1: r[b]=round
        for b in boundary: done.append(b)
        boundary = []
        for n in next: 
            if n in r: boundary.append(n)
#        next=[]

#    for n in G.nodes():
#        if not r.has_key(n): r[n]=round
    return r


def distance_to_nearest_leaf(G,shave_round,shave_limit,trim_degree):
    leaves=[]
    r={}
    for n in G.nodes():
        #print "#",n, G.degree(n),G.degree(n)==1
        if shave_round.get(n,0)<=shave_limit :
            r[n]=0
    for n in G.nodes():
        #print "#",n, G.degree(n),G.degree(n)==1
        if trim_degree.get(n,0)==1 and shave_round.get(n,0)>shave_limit:
            leaves.append(n)
#    print leaves

    round=1
#    print "#",leaves
    boundary=list(leaves)
    next=[]
    done=[]
#    r={}
    for b in boundary: r[b]=round
    while len(boundary)>0:
        round +=1
        next=[]
        for n in boundary:
            for nn in G.neighbors(n):
                if (not nn in boundary+done+next):
                    next.append(nn)
#                    r[nn]=round
        for b in next: 
            r[b] = round
#            if len( [nn for nn in G.neighbors(b) if not r.has_key(nn)] )==1: r[b]=round
        for b in boundary: done.append(b)
        boundary = []
        for n in next: 
            if n in r: boundary.append(n)
#        next=[]
#    for n in G.nodes():
#        if not r.has_key(n): r[n]=round
    return r
    

def shave_round(G):
    leaves=[]
    for n in G.nodes():
        #print "#",n, G.degree(n),G.degree(n)==1
        if G.degree(n)==1:
            leaves.append(n)
#    print leaves

    round=1
#    print "#",leaves
    boundary=list(leaves)
    next=[]
    done=[]
    r={}
    for b in boundary: r[b]=round
    while len(boundary)>0:
        round +=1
        next=[]
        for n in boundary:
            for nn in G.neighbors(n):
                if (not nn in boundary+done+next):
                    next.append(nn)
#                    r[nn]=round
        nr={}
        for b in next: 
            if len( [nn for nn in G.neighbors(b) if nn not in r] )==1: nr[b]=round
        r.update(nr)
        for b in boundary: done.append(b)
        boundary = []
        for n in next: 
            if n in r: boundary.append(n)
#        next=[]
    for n in G.nodes():
        if n not in r: r[n]=round
    return r

def log(x):
    sys.stderr.write(x+"\n")


LOCAL_BRIDGE=1
CUT_ME=2

edge_tags={}

def add_tag(s1,s2,tag):
    if s1<s2:
        ot = edge_tags.get((s1,s2),set())
        ot.add(tag)
        edge_tags[s1,s2] = ot
    else:
        ot = edge_tags.get((s2,s1),set())
        ot.add(tag)
        edge_tags[s2,s1] = ot


def get_tags(s1,s2):
    if s1<s2:
        ot = edge_tags.get((s1,s2),set())
        return tuple(ot)
    else:
        ot = edge_tags.get((s2,s1),set())
        return tuple(ot)


edge_color_setting="hair"
def edge_tag_to_style(tags,setting=edge_color_setting):
    if setting == "hair":
        style=""
        if "hair" in tags:
            style= "color=red"
        elif "longHair" in tags:
            style= "color=orange"
        elif "H" in tags:
            style= "color=blue"
        elif "Y" in tags:
            style= "color=goldenrod"
        elif "nearY" in tags:
            style= "color=goldenrod4"
        elif "bigH" in tags:
            style= "color=green"
        if "promisc" in tags:
            style += " style=dashed"
        return style
            
def printdot(g,gg0,c,n,ll,bh,annot,trim_level={},post_trim_degree={},tag="bad",yDist={},leafDist={},edgeTags=edge_tags):
#def printdot(g,c,n,ll,bh,annot,tag="bad"):
#    print ll

    import colorsys

    chromosomes={}
    chr_mins={}
    chr_maxs={}
    if bh:
        sb=[]
        for cc in c:
            bhi = bh.get(cc,[0,0,0,0,0,0])
            sb.append( ( bhi[1],old_div((float(bhi[3])+float(bhi[4])),2.0),bhi[2],cc ) )
            chromosomes[bhi[1]]=1
            if chr_mins.get(bhi[1],5.0e9)>min( float(bhi[3]), float(bhi[4]) ):  chr_mins[bhi[1]]=min( float(bhi[3]), float(bhi[4]) )
            if chr_maxs.get(bhi[1],-1.0) <max( float(bhi[3]), float(bhi[4]) ):  chr_maxs[bhi[1]]=max( float(bhi[3]), float(bhi[4]) )
        sb.sort()

# {'19': 46194830.0, '18': 59221558.0, '8': 96645227.0, '4': 18548230.0, 'X': 102465955.0}
# {}

    print("#",chr_mins)
    print("#",chr_maxs)

    nchrs=len(list(chromosomes.keys()))
    i=0
    chr_hue={}
    for ch in list(chromosomes.keys()):
        chr_hue[ch] = old_div(float(i),nchrs)
        i+=1

    gg = nx.subgraph(g,c)
    f=open("%s-%d.txt" % (tag,n), "wt")
    nn=1
    lab0={}
    for x in c:
        p=""
        d=x
        if x[0]=="-": 
            p="-"
            d=x[1:]
#        lab0[d]=lab0.get(d,nn)
        lab0[d]=lab0.get( d, float(ll.get(d,0)))
        #print x,d,p,lab0[d]
        nn+=1
    lab={}
    node_fill={}
    for x in c:
        p=""
        d=x
        if x[0]=="-": 
            p="-"
            d=x[1:]

        bhi = bh.get(x,False)
        if bhi:
            lab[x]="{:.1f} {}{}\\n{:.2f}-{:.2f}\\n{}".format( old_div(lab0.get(d,0.0),1000), bhi[1],bhi[2],old_div(float(bhi[3]),1.0e6),old_div(float(bhi[4]),1.0e6), x)
#            lab[x]="{:.1f} {}{}\\n{:.2f}-{:.2f}\\n{} {} {} {}".format( lab0.get(d,0.0)/1000, bhi[1],bhi[2],float(bhi[3])/1.0e6,float(bhi[4])/1.0e6,leafDist.get(x,""),yDist.get(x,""), trim_level[x], post_trim_degree.get(x,"") )
            if ( chr_maxs.get(bhi[1], (1.0+float(bhi[3])+float(bhi[4])) ) ) ==0.0: #/2.0)-chr_mins.get(bhi[1],0.0))==0.0: 
                print("wtf?",x,bhi,bhi[1],( chr_maxs.get(bhi[1], (1.0+float(bhi[3])+float(bhi[4])) ) ))
            rgb=(0,0,0)
            try:
                rgb=colorsys.hls_to_rgb( chr_hue[bhi[1]], 0.5, old_div((old_div((float(bhi[3])+float(bhi[4])),2.0) - chr_mins.get(bhi[1],0)),(chr_maxs.get(bhi[1],old_div((1.0+float(bhi[3])+float(bhi[4])),2.0))-chr_mins.get(bhi[1],0.0)))  )
            except Exception as e:
                print(e)
            node_fill[x]= '#%02x%02x%02x' % (255.0*rgb[0], 255.0*rgb[1], 255.0*rgb[2] ) #"#{}{}{}".format()
        else:
            lab[x]="{:.1f}\\n{}".format(old_div(lab0.get(d,0.0),1000),str(x))
            node_fill[x]="white"
#        lab[x]="{} {}".format( trim_level.get(x,"?"), post_trim_degree.get(x,"?") )

    f.write( "graph G {\n")
#    f.write( "node [margin=0 fontcolor=blue fontsize=32 width=0.5 shape=circle style=filled]")
    f.write( "node [margin=0 fontsize=6 shape=box];\n")
    f.write( "edge [ fontsize=6 ];\n")

    for x in list(lab.keys()):
        f.write( "{0} [label=\"{1}\" fillcolor=\"{2}\" style=\"filled\" color=\"{2}\"] ; \n".format(x,lab[x],node_fill[x]) )

    if bh:
        last=False
        lastx=0.0
        lastc=0
        for c in sb:
            if last and c[0]==last and (c[1]-lastx)<1000000:
                last_bhi = bh.get(lastc,False)
                this_bhi = bh.get(c[-1],False)
                blast_label=str(last_bhi) + str(this_bhi)
                if this_bhi and last_bhi :
                    aa = tuple(last_bhi[1:5])
                    bb = tuple(this_bhi[1:5])
                    qd = qdist(aa,bb)
                    blast_label = "{}".format(qd)
                    
                if gg0.has_edge(lastc,c[-1]) and not t.has_edge(lastc,c[-1]):
                    f.write("\t \"{}\" -- \"{}\" [weight=2 style=dotted label=\"{} {}\" fontcolor=red] ;\n".format(lastc,c[-1],blast_label,int(abs(gg0[lastc][c[-1]]['weight']))))
                else:
                    f.write("\t \"{}\" -- \"{}\" [weight=2 style=dotted label=\"{}\" fontcolor=blue] ;\n".format(lastc,c[-1],blast_label))
            last=c[0]
            lastx=c[1]
            lastc=c[-1]

    for e in gg.edges():
#        f.write( "\t\"%s\" -- \"%s\";\n" % (lab[e[0]],lab[e[1]]) ) #,gg[e[0]][e[1]]['weight'])
#        color="black"
#        if annot.get(e,0)&LOCAL_BRIDGE : color="red"
#        if annot.get(e,0)&CUT_ME       : color="yellow"
        f.write( "\t \"%s\" -- \"%s\" [label=\"%d\" weight=1 %s];\n" % ( e[0],e[1],int(abs(gg[e[0]][e[1]]['weight'])),edge_tag_to_style( get_tags(e[0],e[1]) ) ))

    f.write( "}\n")


def independent_path(G,a,b,k,t):
    q=[a]
    l={}
    l[a]=0
    r=[]
    while len(q)>0:
#        print a,b,G[a][b],q,r
        n=q.pop(0)
        r.append(n)
        for nn in G.neighbors(n):
            if (n==a and nn==b) or (n==b and nn==a): continue
            if G[n][nn]['weight']>-t: continue
            if nn==b: return True
#            print q,[l[i] for i in q]
            l[nn] = min(l.get(nn,10000), l[n]+1)
            if (not nn in q+r) and (l[nn]<=k):
                q.append(nn)
    return False
    
def annotate_edges(t,G,node_list):
    an={}
#    nn= len(list(t.edges()))
#    i=0.0
    for a,b in t.edges():

        if not independent_path(G,a,b,4,2): 
            an[a,b]="local_bridge"
            an[b,a]="local_bridge"
            
    return an


def pairs_overlap(x,y):
    a=min(x[0],x[1])
    b=max(x[0],x[1])
    c=min(y[0],y[1])
    d=max(y[0],y[1])
        
    if a<=c and c<=b: return True
    if a<=d and d<=b: return True
    if c<=a and a<=d: return True
    if c<=b and b<=d: return True
    return False

def qdist(x,y):
    if (not x) or (not y): return (-1)
    if x[0]==y[0]:
        x,y,w,z = list(map(int,[x[2],x[3],y[2],y[3]]))
        ol = pairs_overlap((x,y),(w,z))
        if ol:
            return(-1*min(
                abs(x-w),
                abs(x-z),
                abs(y-w),
                abs(y-z )))
        else:
            return(min(
                abs(x-w),
                abs(x-z),
                abs(y-w),
                abs(y-z )))
    else:
        return(1.0e12)



def pledge_singletons2(g,sg,thresh=2):
    singletons=[]
    ccn=0
    comp={}
    for c in nx.connected_components(sg):
        ccn+=1
        if len(c)==1: singletons.append(c[0])
        for cc in c: comp[cc]=ccn
    for s in singletons:
        total_weight_by_comp={}
        links_by_comp=[]
        exemplar_by_comp={}
        for n in g.neighbors(s):
            ncomp = comp.get(n,-1)
            w = abs(g[s][n]['weight'])
            if w >= thresh:
                links_by_comp.append( (w,ncomp,n) )
                exemplar_by_comp[ncomp]=n

        links_by_comp.sort(reverse=True)
        print("#pledge_stat:",links_by_comp,s)
        if len(links_by_comp)==0: continue
        if len(links_by_comp)==1 or links_by_comp[0][1]==links_by_comp[1][1]:
            n = exemplar_by_comp[links_by_comp[0][1]]
            sg.add_edge( s, n, {'weight': -1} )
            


def pledge_singletons(g,sg,min_combined_weight,min_delta):
    singletons=[]
    ccn=0
    comp={}
    for c in nx.connected_components(sg):
        ccn+=1
        if len(c)==1: singletons.append(c[0])
        for cc in c: comp[cc]=ccn
    for s in sorted(singletons):
        total_weight_by_comp={}
        exemplar_by_comp={}
        for n in sorted(g.neighbors(s)):
            ncomp = comp.get(n,-1)
            total_weight_by_comp[ncomp] = total_weight_by_comp.get(ncomp,0.0) + abs(g[s][n]['weight']  )
            exemplar_by_comp[ncomp]=n
        ncomps = list(total_weight_by_comp.keys())
        ncomps.sort(key=lambda x: total_weight_by_comp[x],reverse=True)
        best = total_weight_by_comp[ncomps[0]]
        delta = 10000.0
        if len(ncomps)>1:
            delta = total_weight_by_comp[ncomps[0]]-total_weight_by_comp[ncomps[1]]
        print("#pledge_stat:",best,delta,s,exemplar_by_comp[ncomps[0]])
        if best >= min_combined_weight and delta >= min_delta:
            sg.add_edge( s, exemplar_by_comp[ncomps[0]], {'weight': -1} )



if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-t','--threshold',default=0.0 ,  type=float)
    parser.add_argument('--define_components_only',default=False,action='store_true')
    parser.add_argument('-H','--head',default=False,type=int)
    parser.add_argument('-D','--outdir',default="link_chunks")
    parser.add_argument('-d','--debug',default=False,action='store_true')

    parser.add_argument('-p','--progress',default=False,action='store_true')
    parser.add_argument('-M','--maxdegree',default=False,type=int)
    parser.add_argument('-m','--minlength',default=500,type=int)

    parser.add_argument('-S','--silent',default=False,action='store_true')
    parser.add_argument('-K','--cutPromisc',default=False,action='store_true')
    parser.add_argument('-J','--logH',default=False,action='store_true')
    parser.add_argument('-T','--logTags',default=False,action='store_true')
    parser.add_argument('-C','--cheat',default=False,action='store_true')
    parser.add_argument('-B','--blacklist')

    parser.add_argument('-b','--besthits')

    parser.add_argument('-c','--nchunks',default=32,type=int)
    parser.add_argument('-l','--lengths')
    parser.add_argument('-E','--edgefile')
    parser.add_argument('-e','--pledgeingedgefile')
    parser.add_argument('-L','--maxLength',type=float,default=150000.0)
    parser.add_argument('-P','--promisc',type=float,default=0.023)
    parser.add_argument('-o','--dotLabel',default="bad")
    parser.add_argument('--seed',required=False,type=int,default=1, help="Seed for random number generation, use -1 for no seed")
    args = parser.parse_args()
    if args.debug:
        args.progress=True
    if args.seed != -1 :
      random.seed(args.seed)
    if args.progress: log( str(args) )
    print("#"+str(args))

    G=nx.Graph()
    SG=nx.Graph()

    ll={}
    if args.lengths:
        f = open(args.lengths)
        while True:
            l = f.readline()
            if not l: break
            if l[0]=="#": continue

            c=l.strip().split()
            l = int(c[1])
            ll[c[0]]=int(c[1])
            if l>= args.minlength:
                G.add_node(c[0])
                SG.add_node(c[0])
        f.close()
    if args.progress: print("#Done reading lengths")



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


    if args.edgefile:
        f = open(args.edgefile)
    else:
        f=sys.stdin
    while True:
        l = f.readline()
        if not l: break
        if l[0]=="#": continue
        c=l.strip().split()
        u,v,w = c[0],c[1],float(c[2])
        if ( not args.lengths ) or (ll[u]>=args.minlength and ll[v]>=args.minlength):
            G.add_edge(u,v,weight=-w)
            SG.add_node(u)
            SG.add_node(v)
            if w >= args.threshold:
                SG.add_edge(u,v,weight=-w)
    if args.edgefile:
        f.close()

    if args.pledgeingedgefile:
        f = open(args.pledgeingedgefile)
        while True:
            l = f.readline()
            if not l: break
            if l[0]=="#": continue
            c=l.strip().split()
            u,v,w = c[0],c[1],float(c[2])
            if ( not args.lengths ) or (ll[u]>=args.minlength and ll[v]>=args.minlength):
                G.add_edge(u,v,weight=-w)
        f.close()
        
    if args.progress: print("#Done reading edgelist")

    bad_nodes=[]
    total_discarded_length=0
    total_discarded_length1=0
    n_discarded1=0
    n_discarded2=0
    if args.maxdegree:
        for n in sorted(SG.nodes()):
            print("#dg:", SG.degree(n))
            if SG.degree(n)>args.maxdegree:
                n_discarded1+=1
                bad_nodes.append(n)
                total_discarded_length += ll[n]
                print("#discard:",n,ll[n],SG.degree(n))
                for nn in SG.neighbors(n):
                    if SG.degree(nn)==1:
                        n_discarded2+=1
                        total_discarded_length1+=ll[nn]
        for n in bad_nodes:
            e_to_remove=[]
            for e in SG.edges([n]):
                e_to_remove.append(e)
            SG.remove_edges_from(e_to_remove)

    if args.cutPromisc:
        bad_nodes=[]
        for n in sorted(SG.nodes()):
            print("#ps:", old_div(float(SG.degree(n)), ll[n]), args.promisc,SG.degree(n),ll[n])
            print("#pr:", old_div(float(G.degree(n)), ll[n]), args.promisc,G.degree(n),ll[n])
            if (old_div(float(G.degree(n)), ll[n]))>args.promisc : # G.degree(n)/ll[n]>args.maxdegree:
                n_discarded1+=1
                bad_nodes.append(n)
                total_discarded_length += ll[n]
                print("#discard:",n,ll[n],G.degree(n))
#                for nn in G.neighbors(n):
#                    if G.degree(nn)==1:
#                        n_discarded2+=1
#                        total_discarded_length1+=ll[nn]
        for n in bad_nodes:
            e_to_remove=[]
            for e in SG.edges([n]):
                e_to_remove.append(e)
            #G.remove_edges_from(e_to_remove)
            SG.remove_edges_from(e_to_remove)
             
    if args.blacklist:
        f=open(args.blacklist)
        e_to_remove=[]
        while True:
            l=f.readline()
            if not l: break
            c=l.strip().split()
            e_to_remove.append((c[1],c[2]))
        G.remove_edges_from(e_to_remove)
        SG.remove_edges_from(e_to_remove)
        
            
    if args.cheat:
        e_to_remove=[]
        for a,b in sorted(SG.edges()):
            if not ( a in besthit and b in besthit):
                e_to_remove.append((a,b))
            else:
                aa = tuple(besthit[a][1:5])
                bb = tuple(besthit[b][1:5])
                qd = qdist(aa,bb)
                if qd >= args.maxLength : 
                    e_to_remove.append((a,b))
        SG.remove_edges_from(e_to_remove)
        

    if args.progress: print("#total_discarded_length",n_discarded1,n_discarded2,old_div(float(total_discarded_length),1.0e6),old_div(float(total_discarded_length1),1.0e6),old_div(float(total_discarded_length+total_discarded_length1),1.0e6))


    promisc = {}
    for n in sorted(SG.nodes()):
        if args.debug: print("#r:",old_div(float(SG.degree(n)),ll[n]))
        if (old_div(float(SG.degree(n)),ll[n]))>args.promisc: promisc[n]=True

    tag_tallies={}
    bad_tag_tallies={}
    strx={"+":0, "-":1}
    strings = []
    ccn=0
    component={}

    component_contigs={}
    chunk={}

    #pledge_singletons(G,SG,args.threshold,2)
    pledge_singletons2(G,SG,2)


    for c in sorted(nx.connected_components(SG), key=lambda x: " ".join(sorted(x))) :
        ccn+=1
        print("c:",ccn,ccn%args.nchunks,len(c),sorted(c))
        for cc in c: 
            component[cc]=ccn
            chunk[cc]=ccn%args.nchunks
        component_contigs[ccn] = tuple(c)

    if args.define_components_only: exit(0)

    intra_fhs={}
    inter_fhs={}
    for i in range(args.nchunks):
        intra_fhs[i]       = open("{}/intra.{}.links".format(args.outdir,i),"wt")
        for j in range(i,args.nchunks):
            inter_fhs[i,j] = open("{}/inter.{}-{}.links".format(args.outdir,i,j),"wt")
            inter_fhs[j,i] = inter_fhs[i,j]            

    while True:
        l = sys.stdin.readline()
        if not l: break
        if l[0]=="#": continue
        c=l.strip().split("\t")
        if not (ll[c[0]]>=args.minlength and ll[c[1]]>=args.minlength):  continue
#        if ( not args.lengths ) or (ll[u]>=args.minlength and ll[v]>=args.minlength):

        if args.debug: 
            try:
                print("#",c[0],c[1],component.get(c[0]),component.get(c[1]))
            except Exception as e:
                print(e)
                print("#wtf",l)
        chunk1,chunk2 = chunk[c[0]],chunk[c[1]]
        comp1,comp2 = component[c[0]],component[c[1]]
        if comp1==comp2:
            intra_fhs[chunk1].write(l)
        else:
            inter_fhs[chunk1,chunk2].write(l)
            
#        if component.get(c[0]) and component.get(c[1]) and component.get(c[0])==component.get(c[1]):
#            print l.strip()
        

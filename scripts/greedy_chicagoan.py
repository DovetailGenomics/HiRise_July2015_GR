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

default_gapsize=100.0

ll={}
links={}

def is_shaved_tail(G,shave_round,shaved_degree,shave_limit):
    leaves=[]
    r={}
    for n in G.nodes():
        if shave_round.get(n,0)>shave_limit and shaved_degree.get(n,0)==1:
            leaves.append(n)
    for l in leaves:
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

def test_inversion_option(c1,c2,join_options,ograph,linked):
    a=c1
    b=linked[c1]
    c=c2
    d=linked[c2]

    #k=linked[internal_node]
    ograph.remove_edge(a,b)
    ograph.remove_edge(c,d)
                                                    ##                                                                        ---b a-------c d----
    if a in nx.node_connected_component(ograph,c):  ## exchance labels of a,b if necessary, so the nodes are in this config:  ---a b-------c d------ 
        x=b
        b=a
        a=x
     
                                                    ##                                                                        ---b a-------d c----
    if a in nx.node_connected_component(ograph,d):  ## exchance labels of a,b if necessary, so the nodes are in this config:  ---a b-------c d------ 
        x=b
        b=a
        a=x

        x=d
        d=c
        c=x

                                                    ##                                                                        ---a b-------d c----
    if b in nx.node_connected_component(ograph,d):  ## exchance labels of c,d if necessary, so the nodes are in this config:  ---a b-------c d------ 
        x=d
        d=c
        c=x

    n_scaffold = old_div(len(nx.node_connected_component(ograph,b)),2)
    print("inversion n nodes",n_scaffold)


    total_i_len = sum( ograph[i][j]['length'] for i,j in nx.bfs_edges(ograph,b) )
    print("inv len",total_i_len)
    if total_i_len < 10000.0 or n_scaffold<2:
        print("inversion length",total_i_len,n_scaffold,"too short")
        join_options.append( (0.0,(),() ) )                                                                                              
        ograph.add_edge(a,b,length=default_gapsize,contig=False)
        ograph.add_edge(c,d,length=default_gapsize,contig=False)

        return

    interc_score0 = intercalation_score_raw(a,b,d,ograph) 
    interc_score1 = intercalation_score_raw(a,c,d,ograph) 
    print("inversion0", interc_score0)
    print("inversion", interc_score1)
    join_options.append( (interc_score0,(),() ) )
    join_options.append( (interc_score1,((a,c),(b,d)),((a,b),(c,d)) ) )
    ograph.add_edge(a,b,length=default_gapsize,contig=False)
    ograph.add_edge(c,d,length=default_gapsize,contig=False)

    return

def test_endInversion_option(free_end,internal_node,join_options,ograph,linked):
    if free_end in linked:
        x=free_end
        free_end = internal_node
        internal_node = x
    print("end inversion",free_end,linked.get(free_end),internal_node,linked[internal_node])

    k=linked[internal_node]
    ograph.remove_edge(internal_node,k)

    if free_end in nx.node_connected_component(ograph,internal_node): 
        x=k
        k=internal_node
        internal_node=x

    sc = link_test(ograph,k,internal_node)
    join_options.append( (sc,(),()))
    print("end inversion existing:",sc)

    sc = link_test(ograph,free_end,internal_node)
    join_options.append( (sc ,((free_end,internal_node),),((internal_node,k),) ) ) 
    print("end inversion:",sc)

    ograph.add_edge(internal_node,k,length=default_gapsize, contig=False) 
    return

def test_interc_option( gap_edge1, gap_edge2, free_end1, join_options, ograph ):
    if not ograph.has_edge(gap_edge1,gap_edge2):
        print("expected nodes to be connected: {} {}".format(gap_edge1, gap_edge2))
        raise Exception('not connected')
    ograph.remove_edge( gap_edge1, gap_edge2 )
    if gap_edge1 in nx.node_connected_component(ograph,gap_edge2): 
        print("problem: these should be disconnected now:",gap_edge1,gap_edge2)
        raise Exception('too connected i')
    if gap_edge2 in nx.node_connected_component(ograph,free_end1): 
        print("problem: these should have been disconnected all along",free_end1,gap_edge2)
        raise Exception('too connected j')
    if gap_edge1 in nx.node_connected_component(ograph,free_end1): 
        print("problem: these should have been disconnected all along",free_end1,gap_edge1) 
        raise Exception('too connected k')
    interc_score = intercalation_score(gap_edge1,free_end1,gap_edge2,ograph)
    ograph.add_edge(gap_edge1, gap_edge2, length=default_gapsize,contig=False)
    d = far_end(ograph,free_end1)
    join_options.append( (interc_score,((gap_edge1,free_end1),(d,gap_edge2)),((gap_edge1, gap_edge2),) ) )


N=76710553.0 #100000000.0
pn=0.3
GenomeSize=3.0e9
import math

def get_score(c1_stranded,c2_stranded,gaplen=default_gapsize,cache={}):
    if gaplen<0.0: print(math.log(-1))
#    if cache.has_key((c1_stranded,c2_stranded,gaplen)):
#        return cache[c1_stranded,c2_stranded,gaplen]
#    if cache.has_key((c2_stranded,c1_stranded,gaplen)):
#        return cache[c2_stranded,c1_stranded,gaplen]
    o1,o2=0,0
    if c1_stranded[-2:]==".5": o1=1
    if c2_stranded[-2:]==".3": o2=1
    c1=c1_stranded[:-2]
    c2=c2_stranded[:-2]
    l1=ll[c1]
    l2=ll[c2]
    p0 = ces.p_not_a_hit(l1,l2,GenomeSize,gaplen,pn)
    if p0<0.0:
        print("wtf?  p0<0",p0,l1,l2,GenomeSize,gaplen,pn)
    thisscore=ces.llr_v0(l1,l2,o1,o2,GenomeSize,pn,links.get((c1,c2),[]),N,gaplen,p0 )
    #cache[c1_stranded,c2_stranded,gaplen]=thisscore
    return thisscore

def traverse_and_layout(n,coords,facing,x,s,og,maxD=2000000):
    color={}
    q=[n]
    ptr=x
    if s>0: 
        facing[n]="L"    
    else: 
        facing[n]="R"
    while (abs(ptr-x)<maxD) and len(q)>0:
        m=q.pop(0)
        coords[m]=ptr
        color[m]=1
        for mm in og.neighbors(m):
            if not mm in color:
                q.append(mm)
                ptr+=s*og[m][mm]['length']
                if not og[m][mm]['contig']:
                   if s==-1:
                       facing[mm]='R'
                   else:
                       facing[mm]='L'
                else:
                   if s==-1:
                       facing[mm]='L'
                   else:
                       facing[mm]='R'

def far_end(g,n):
    if not g.degree(n)==1:
        print("wtf: this should be a leaf")
        exit(0)
    for m in nx.node_connected_component(g,n):
        if (not m==n) and g.degree(m)==1:
            return m

def both_ends(g,n):
    p=[]
    for m in nx.node_connected_component(g,n):
        if  g.degree(m)==1:
            p.append( m)
    return p

def strand_check(a,b,x,f):
    if x[a]<x[b]:
        if f[a]=="R" and f[b]=="L": return True
    else:
        if f[b]=="R" and f[a]=="R": return True
    return False

def take_best_join_option(join_options,og,linked,edge_options,label=""):
#    print join_options
    join_options.sort(reverse=True)
#    print join_options
    if join_options[0][0]>0.0:
#        if len(join_options[0][2]) >=1 : print "#accept join",label, join_options[0], join_options
        for x,y in join_options[0][2]:
            og.remove_edge(x,y)
            del linked[x]
            del linked[y]
        for x,y in join_options[0][1]:
            linked[x]=y
            linked[y]=x
            og.add_edge(x,y,length=default_gapsize,contig=False)
            edge_options[x,y] = "[color=brown label=\"{:.2f}\"]".format(weights.get((x,y),0))
            edge_options[y,x] = "[color=brown label=\"{:.2f}\"]".format(weights.get((x,y),0))
                    
def link_test(og,a,b,gapsize=default_gapsize,max_stretch=200000):
    coords = {}
    facing={}
    traverse_and_layout(a,coords,facing,0,-1,og,maxD=200000)
    traverse_and_layout(b,coords,facing,gapsize,+1,og,maxD=200000)
#    print "#x:setup done"
    sys.stdout.flush()
    
    score=0.0
    for n1 in nx.node_connected_component(og,a):
        if abs(coords.get(n1,-1e10))>max_stretch: continue
        for n2 in nx.node_connected_component(og,b):
            if abs(coords.get(n2,1e10))>max_stretch: continue
            if strand_check(n1,n2,coords,facing):
                distance = coords[n2] - coords[n1]
                if distance < max_stretch:
                    x=get_score(n1,n2,distance)
#                    print "#x:partial",n1,n2,x
                    sys.stdout.flush()
                    score += x
    return score

def intercalation_score(a,b,c,og):
#    print "#intercalation test:",a,b,c
    if a in nx.node_connected_component(og,b): 
        print("a should not be connected to b",a,b)
        raise Exception('too connected a')
    if b in nx.node_connected_component(og,c): 
        print("b should not be connected to c",b,c)
        raise Exception('too connected b')
    if a in nx.node_connected_component(og,c): 
        print("a should not be connected to c",a,c)
        raise Exception('too connected c')
    coordinates1={}
    coordinates2={}
    facing={}
    traverse_and_layout(a,coordinates1,facing,0,-1,og)
    coordinates2=dict(coordinates1)
    traverse_and_layout(c,coordinates1,facing,default_gapsize,+1,og)
    traverse_and_layout(b,coordinates2,facing,default_gapsize,+1,og)
    traverse_and_layout(c,coordinates2,facing,default_gapsize+max(coordinates2.values()),+1,og)

    score_0=0.0
    for n1 in nx.node_connected_component(og,a):
        for n2 in nx.node_connected_component(og,c):
            if strand_check(n1,n2,coordinates1,facing):
                distance = coordinates1[n2] - coordinates1[n1]
                if distance < 200000:
                    if distance<0: print("wtf1?")
                    score_0 += get_score(n1,n2,distance)

    score_1=0.0
    for n1 in nx.node_connected_component(og,a):
        for n2 in nx.node_connected_component(og,c):
            if strand_check(n1,n2,coordinates2,facing):
                distance = coordinates2[n2] - coordinates2[n1]
                if distance < 200000:
                    if distance<0: print("wtf2?")
                    score_1 += get_score(n1,n2,distance)

        for n2 in nx.node_connected_component(og,b):
            if strand_check(n1,n2,coordinates2,facing):
                distance = coordinates2[n2] - coordinates2[n1]
                if distance < 200000:
                    if distance<0: print("wtf3?")
                    score_1 += get_score(n1,n2,distance)

    for n1 in nx.node_connected_component(og,b):
        for n2 in nx.node_connected_component(og,c):
            if strand_check(n1,n2,coordinates2,facing):
                distance = coordinates2[n2] - coordinates2[n1]
                if distance < 200000:
                    if distance<0: print("wtf4?")
                    score_1 += get_score(n1,n2,distance)

#    print "#intercalation scores:",score_1,score_0


    if a in nx.node_connected_component(og,b): 
        print("a should not be connected to b",a,b)
        raise Exception('too connected d')
    if b in nx.node_connected_component(og,c): 
        print("b should not be connected to c",b,c)
        raise Exception('too connected e')
    if a in nx.node_connected_component(og,c): 
        print("a should not be connected to c",a,c)
        raise Exception('too connected f')

    return (score_1-score_0)

def intercalation_score_raw(a,b,c,og):
#    print "#intercalation test:",a,b,c
    if a in nx.node_connected_component(og,b): print("a should not be connected to b",a,b)
    if b in nx.node_connected_component(og,c): print("b should not be connected to c",b,c)
    if a in nx.node_connected_component(og,c): print("a should not be connected to c",a,c)
#    coordinates1={}
    coordinates2={}
    facing={}
    traverse_and_layout(a,coordinates2,facing,0,-1,og)
#    coordinates2=dict(coordinates1)
#    traverse_and_layout(c,coordinates1,facing,1000,+1,og)
    traverse_and_layout(b,coordinates2,facing,default_gapsize,+1,og)
    traverse_and_layout(c,coordinates2,facing,default_gapsize+max(coordinates2.values()),+1,og)

#    score_0=0.0
#    for n1 in nx.node_connected_component(og,a):
#        for n2 in nx.node_connected_component(og,c):
#            if strand_check(n1,n2,coordinates1,facing):
#                distance = coordinates1[n2] - coordinates1[n1]
#                if distance < 200000:
#                    if distance<0: print "wtf1?"
#                    score_0 += get_score(n1,n2,distance)

    score_1=0.0
    for n1 in nx.node_connected_component(og,a):
        for n2 in nx.node_connected_component(og,c):
            if strand_check(n1,n2,coordinates2,facing):
                distance = coordinates2[n2] - coordinates2[n1]
                if distance < 200000:
                    if distance<0: print("wtf2?")
                    score_1 += get_score(n1,n2,distance)

        for n2 in nx.node_connected_component(og,b):
            if strand_check(n1,n2,coordinates2,facing):
                distance = coordinates2[n2] - coordinates2[n1]
                if distance < 200000:
                    if distance<0: print("wtf3?")
                    score_1 += get_score(n1,n2,distance)

    for n1 in nx.node_connected_component(og,b):
        for n2 in nx.node_connected_component(og,c):
            if strand_check(n1,n2,coordinates2,facing):
                distance = coordinates2[n2] - coordinates2[n1]
                if distance < 200000:
                    if distance<0: print("wtf4?")
                    score_1 += get_score(n1,n2,distance)

#    print "#intercalation scores:",score_1,score_0
    return (score_1)


if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()
   
    parser.add_argument('-t','--threshold',default=0.0 ,  type=float)

    parser.add_argument('-H','--head',default=False,type=int)
    parser.add_argument('-D','--savetreedots',default=False,action='store_true')
    parser.add_argument('-d','--debug',default=False,action='store_true')
#    parser.add_argument('-I','--nointerc',default=False,action='store_true')
    parser.add_argument('-p','--progress',default=False,action='store_true')
    parser.add_argument('-M','--maxdegree',default=False,type=int)
    parser.add_argument('-m','--minlength',default=500,type=int)

    parser.add_argument('-S','--silent',default=False,action='store_true')
    parser.add_argument('-K','--cutPromisc',default=False,action='store_true')
    parser.add_argument('-J','--logH',default=False,action='store_true')
    parser.add_argument('-T','--logTags',default=False,action='store_true')
    parser.add_argument('-C','--cheat',default=False,action='store_true')
    parser.add_argument('-B','--blacklist')

    parser.add_argument('-j','--joins') #pre-join these
    parser.add_argument('-k','--links') 
    parser.add_argument('-b','--besthits')
    parser.add_argument('-s','--skip',default=0,type=int) 
    parser.add_argument('-l','--lengths')
    parser.add_argument('-E','--edgefile')
    parser.add_argument('--set_insert_size_dist_fit_params')
    parser.add_argument('-L','--maxLength',type=float,default=150000.0)
    parser.add_argument('-P','--promisc',type=float,default=0.023)
    parser.add_argument('-o','--dotLabel',default="bad")

    args = parser.parse_args()
    if args.debug:
        args.progress=True

    if args.progress: log( str(args) )
    print("#"+str(args))

    if args.set_insert_size_dist_fit_params:
        s=args.set_insert_size_dist_fit_params
        a,b,c,d,f,pn,N = list(map(float,list(s.split(','))))
        print("#",a,b,c,d,f)
        ces.set_insert_size_dist_fit_params(a,b,c,d,f,pn,N)
        sys.stdout.flush()

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

    #broken.al.masked.links.txt
    def read_links(contigs=False):
        links={}
        import glob
        for f in glob.glob(args.links): #["same_component-t4.links.txt"]:
            print("#file:",f)
            ff=open(f)
            while True:
                l = ff.readline()
                if not l: break
                if l[0]=="#": continue
                c=l.strip().split("\t")
                if (not contigs) or (c[0] in contigs) and (c[1] in contigs):
                    links[(c[0],c[1])]=eval(c[5])
        return links


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
            if w >= args.threshold:
                SG.add_edge(u,v,weight=-w)
    if args.edgefile:
        f.close()

    if args.progress: print("#Done reading edgelist")

    bad_nodes=[]
    total_discarded_length=0
    total_discarded_length1=0
    n_discarded1=0
    n_discarded2=0
    if args.maxdegree:
        for n in SG.nodes():
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
        for n in G.nodes():
#            print "#dg:", SG.degree(n)
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
            for e in G.edges([n]):
                e_to_remove.append(e)
            G.remove_edges_from(e_to_remove)
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
        for a,b in SG.edges():
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
    for n in G.nodes():
        if args.debug: print("#r:",old_div(float(G.degree(n)),ll[n]))
        if (old_div(float(G.degree(n)),ll[n]))>args.promisc: promisc[n]=True

    tag_tallies={}
    bad_tag_tallies={}
    strx={"+":0, "-":1}
    strings = []
    ccn=1

    if args.progress: 
        print("about to load the links")
        sys.stdout.flush()
    links = read_links()
    if args.progress: 
        print("loaded the links")
        sys.stdout.flush()


    print("#n_connected_components:",len(nx.connected_components(SG)))
    for c in nx.connected_components(SG):
        if ccn<=args.skip:
            ccn+=1
            continue

        if len(c)==1: 
            print("#edge:",c[0]+".5",c[0]+".3",{"length": ll[c[0]],"contig": True})
            continue

#        if args.progress: 
#            print "about to load the links"
#            sys.stdout.flush()
#        links = read_links(c)
#        if args.progress: 
#            print "loaded the links"
#            sys.stdout.flush()

        print("cc:",ccn,len(c),c)

        edge_options={}
        nodes={}

        og=nx.Graph()
        linked={}

        for cc in c:
            og.add_node(cc+".5")
            og.add_node(cc+".3")
            nodes[cc+".5"]="5'"
            nodes[cc+".3"]="3'"
            og.add_edge(cc+".5",cc+".3",length=ll[cc],contig=True)
            edge_options[cc+".5",cc+".3"] = "[label=\"{} {}\" penwidth=3]".format(cc,old_div((ll[cc]+500),1000))
            edge_options[cc+".3",cc+".5"] = "[label=\"{} {}\" penwidth=3]".format(cc,old_div((ll[cc]+500),1000))

        if args.joins:
            f = open(args.joins)
            while True:
                l = f.readline()
                if not l: break
                cz=l.strip().split()
                if not len(cz)>=3: continue
                #print "z:",cz[1][:-2],cz[2][:-2],cz
                if cz[0]=="#join:" :
                    if (cz[1][:-2] in c) and (cz[2][:-2] in c):
                        print("join",cz[1:])
                        og.add_edge( cz[1],cz[2],length=default_gapsize,contig=False)
                        linked[cz[1]]=cz[2]
                        linked[cz[2]]=cz[1]
#                        og.add_edge(c1,c2,length=1000.0,contig=False)
                        edge_options[cz[1],cz[2]] = "[label=\"pre\" penwidth=2]"
                        edge_options[cz[2],cz[1]] = "[label=\"pre\" penwidth=2]"
                    else: 
                        pass
#                        print "skip",cz[1:]
            f.close()

        weights = {}
        for c1 in c:
            if args.progress: print("#",c1)
            l1=ll[c1]
            for c2 in c:
                if not c1<c2: continue
                if args.joins and c2+".3" in nx.node_connected_component(og, c1+".3"): continue
                l2=ll[c2]
                if (c1,c2) not in links: continue
                for gaplen in [default_gapsize]: # [ 0, 500, 1000, 2000 , 10000, 20000, 30000, 40000, 100000, 200000, 500000 ]:
                    p0 = ces.p_not_a_hit(l1,l2,GenomeSize,gaplen,pn) 
                    s={}
                    for (o1,o2,suf1,suf2) in ((0,0,".3",".5"),(0,1,".3",".3"),(1,0,".5",".5"),(1,1,".5",".3")):
                        weights[c1+suf1,c2+suf2] = ces.llr_v0( l1,l2,o1,o2,GenomeSize,pn,links.get((c1,c2),[]),N,gaplen,p0 )
                        weights[c2+suf2,c1+suf1] = weights[c1+suf1,c2+suf2] 
                        if args.progress: print("#",c1+suf1,c2+suf2,len(links[c1,c2]),weights[c1+suf1,c2+suf2])

        link_pairs = list(weights.keys())
        link_pairs.sort(key = lambda x: weights[x], reverse=True)

        for c1,c2 in link_pairs:
            if not c1<c2: continue
            if weights[c1,c2]>12:
                if c2 in nx.node_connected_component(og, c1):  # link within one of the scaffolds.  
                    if c1 in linked and linked[c1]==c2: continue # this join pre made!
                    if (not c1 in linked) and (not c2 in linked):  # this would circularize
                        pass
                    elif ( c1 in linked) and ( c2 in linked):      # test to invert (or excize circle?)
                        if og.has_edge(c1,c2): 
                            print("skip test flip of",c1,c2)
                            continue #don't test inversions of individual contigs
#                        print "test invert",c1,linked.get(c1),c2, linked.get(c2) 
                        join_options=[]
                        test_inversion_option(c1,c2,join_options,og,linked)
                        take_best_join_option(join_options,og,linked,edge_options,"inversion")
                        
                    elif (not c1 in linked) or (not c2 in linked): # test for end inversion (scorpion tail?) or pinch off circle from end?
#                        print "test end invert",c1,c2,linked.get(c1), linked.get(c2)
                        join_options=[]
                        test_endInversion_option(c1,c2,join_options,og,linked)
                        take_best_join_option(join_options,og,linked,edge_options,"end inversion")                        
                        
                else:
                    if (not c1 in linked) and (not c2 in linked): # and not c2 in nx.node_connected_component(og, c1):
                        linked[c1]=c2
                        linked[c2]=c1
                        print("#easy",c1,c2,weights[c1,c2])
                        og.add_edge(c1,c2,length=default_gapsize,contig=False)
                        edge_options[c1,c2] = "[color=blue label=\"{:.2f}\"]".format(weights[c1,c2])
                        edge_options[c2,c1] = "[color=blue label=\"{:.2f}\"]".format(weights[c1,c2])
                    elif (c1 in linked) and (c2 in linked):
                        # This edge looks like an H.  test whether the two linked scaffolds can be put end-to-end:
                        # a-b k-d : we'll test b+k, a+k, b+d and a+d, and if the best one raises the total score, we'll add it.
                        a,b = both_ends(og,c1)
                        k,d = both_ends(og,c2)
                        sc={}
                        for x,y in ((b,k),(a,k),(b,d),(a,d)):
                            sc[x,y] = link_test(og,x,y)
                        sll=list(sc.items())
                        sll.sort(key=lambda x: x[1],reverse=True)
                        print("# Htest:", c1,c2,sll)
                        x,y = sll[0][0]
                        if sc[x,y]>0:
                            print("#H accept", sc[x,y] ,x,y)
                            linked[x]=y
                            linked[y]=x
                            og.add_edge(x,y,length=default_gapsize,contig=False)
                            edge_options[x,y] = "[color=green label=\"{:.2f}\"]".format(weights.get((x,y),0))
                            edge_options[y,x] = "[color=green label=\"{:.2f}\"]".format(weights.get((x,y),0))
                            
                    elif (not c1 in linked) or (not c2 in linked):  #test whether one can be intercalated in the other:
                        accepted=False
                        join_options=[]
                        if (not c1 in linked):
                            b=c1
                            a=c2
                            cP=linked[a]
                            print("#1:",c1,c2,a,b,cP,weights[c1,c2])
                        elif (not c2 in linked):
                            b=c2
                            a=c1
                            cP=linked[a]
                            print("#2:",c1,c2,a,b,cP,weights[c1,c2])
                        if (not c1 in linked) or (not c2 in linked):
                            # not linking end-to-end:  
                            #but one of them IS an end link, try inserting that scaffold into the other one.
                            # a and c are the contig ends flanking the gap we're inserting into, b is the end of
                           #  the scaffold being tested for insertion.


                            test_interc_option(a,cP,b,join_options,og)
                            test_interc_option(a,cP,far_end(og,b),join_options,og)

                            a,k = both_ends(og,a)
                            b,d = both_ends(og,b)
                            sc={}
                            for x,y in ((a,b),(k,b),(a,d),(k,d)):
                                sc[x,y] = link_test(og,x,y)
                                join_options.append( (sc[x,y],((x,y),),() ) )

                            take_best_join_option(join_options,og,linked,edge_options,"interc or end")
                                    

        if False:
            for c1,c2 in link_pairs:
                if weights[c1,c2]>12 and not c2 in nx.node_connected_component(og, c1):
    #                linked[c1]=True
    #                linked[c2]=True
                    og.add_edge(c1,c2)        
                    edge_options[c1,c2] = "[style=dotted label=\"{:.2f}\"]".format(weights[c1,c2]) #]"
                    edge_options[c2,c1] = "[style=dotted label=\"{:.2f}\"]".format(weights[c1,c2]) #]"

        red_edges=""
        node_fill={}
        if args.besthits:
            import colorsys

            chromosomes={}
            chr_mins={}
            chr_maxs={}
            ends=[]
            
            for cc in c:
#best:   Scaffold91_1    46380.0 chr12   +       49064574        49095652    
#                besthit[c[1]]=c[2:]

                if cc not in besthit: continue
                aa = tuple(besthit[cc][1:5])
                chromosomes[aa[0]]=1
                xxx = old_div((int(aa[2])+int(aa[3])),2)
                if chr_mins.get(aa[0],1.0e9)>xxx: chr_mins[aa[0]]=xxx
                if chr_maxs.get(aa[0],-1.0e9)<xxx: chr_maxs[aa[0]]=xxx
                if aa[1]=="+":
                    ends.append((aa[0],xxx,cc+".5",cc+".3",cc))
#                    ends.append((c[0],int(aa[3]),cc+".3",cc))
                elif aa[1]=="-":
                    ends.append((aa[0],xxx,cc+".3",cc+".5",cc))
#                    ends.append((c[0],int(aa[3]),cc+".5",cc))
            ends.sort()

            nchrs=len(list(chromosomes.keys()))
            if nchrs==0: nchrs+=1
            i=0
            chr_hue={}
            for ch in list(chromosomes.keys()):
                chr_hue[ch] = old_div(float(i),nchrs)
                i+=1

            for cc in c:
                if not cc in besthit: continue
                aa = tuple(besthit[cc][1:5])
                chrid=aa[0]
                xxx = old_div(float(int(aa[2])+int(aa[3])),2)
                rgb=colorsys.hls_to_rgb( chr_hue[chrid], 0.5, old_div((xxx - chr_mins.get(chrid,0)),(1.0 + chr_maxs.get(chrid,xxx+1.0)-chr_mins.get(chrid,0.0))))
                node_fill[cc+".5"]= '#%02x%02x%02x' % (255.0*rgb[0], 255.0*rgb[1], 255.0*rgb[2] ) #"#{}{}{}".format()
                node_fill[cc+".3"]= '#%02x%02x%02x' % (255.0*rgb[0], 255.0*rgb[1], 255.0*rgb[2] ) #"#{}{}{}".format()



            for i in range(1,len(ends)):
                if ends[i-1][0]==ends[i][0] and ends[i-1][3] in nx.node_connected_component(og,ends[i][3]):
#                    red_edges += "\"{}\" -- \"{}\" [color=red style=dotted label=\"{:.1f}\" fontcolor=red];\n".format(ends[i-1][3],ends[i][2],weights.get((ends[i-1][3],ends[i][2]),""))
                    red_edges += "\"{}\" -- \"{}\" [color=red style=dotted ];\n".format(ends[i-1][3],ends[i][2])

        if args.savetreedots:
            f=open("greedy-linear7-{}.dot".format(ccn),"wt")
            f.write("graph G {\n")
            for cc in list(nodes.keys()):
    #[label=\"{1}\" fillcolor=\"{2}\" style=\"filled\" color=\"{2}\"]
                f.write("\"{}\" [ label=\"{}\" fillcolor=\"{}\" color=\"{}\" style=\"filled\" ];\n".format(cc,nodes[cc],node_fill.get(cc,"white"),node_fill.get(cc,"white")))
                #            f.write("\n".format(c))
            f.write(red_edges)
            for c1,c2 in og.edges():
                f.write("\"{}\" -- \"{}\" {} \n".format(c1,c2,edge_options.get((c1,c2),"")))
            f.write("}\n")
            f.close()
      
        for nnog in nx.connected_component_subgraphs(og):
            tlnnog=0
            for c1,c2 in nnog.edges():
                print("#edge:",c1,c2,og.get_edge_data(c1,c2))
                tlnnog+=og[c1][c2]['length']
            print("#slen",tlnnog)
        ccn+=1

        if args.head and ccn>args.head:
            break

    



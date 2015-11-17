#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import map
from builtins import range
from past.utils import old_div
import sys
import networkx as nx
#import chicago_edge_scores as ces
#import orienterN4 as orient
import local_oo_opt as loo
#import greedy_chicagoan gc

default_gapsize=100.0

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
        if "blacklist" in tags:
            style="penwidth=5 color=\"red\""
            return style

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
            bhi = bh.get(cc,["0","0","0",0,0,"0"])
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

    parser.add_argument('-c','--mychunk',type=int,default=0)
    parser.add_argument('-b','--besthits')
    parser.add_argument('-l','--lengths')
    parser.add_argument('-Z','--chunkfile')
    parser.add_argument('-E','--edgefile')
    parser.add_argument('-k','--links',default=False)
    parser.add_argument('-L','--maxLength',type=float,default=150000.0)
    parser.add_argument('-P','--promisc',type=float,default=0.023)
    parser.add_argument('-o','--dotLabel',default="bad")
    parser.add_argument('--set_insert_size_dist_fit_params')

    args = parser.parse_args()
    if args.debug:
        args.progress=True

    if args.progress: log( str(args) )
    print("#"+str(args))


    if True:
        
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
            loo.set_fit_params(fit_params)
        
    blacklist={}
    if args.blacklist:
        for l in open(args.blacklist):
            if l[0]=="#": continue
            blacklist[tuple(l.strip().split()[:2])]=1
            
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
    #orient.length=ll

    links={}
    if args.links:
        f=open(args.links)
        while True:
            l = f.readline()
            if not l: break
            if l[0]=="#": continue
            c = l.strip().split()
            if len(c)<5: continue
            
            pp = eval(" ".join(c[5:]))
            links[c[0],c[1]]= pp
            links[c[1],c[0]]= [ (p[1],p[0]) for p in pp ] 
            G.add_edge(c[0],c[1],weight=-float(c[4]))
            
        f.close()
    if args.progress: print("#Done reading links", len(G.nodes()), len(G.edges()))

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

    components=[]
    component={}
    
    if args.chunkfile:
        f = open(args.chunkfile)
        while True:
            l = f.readline()
            if not l: break
            c = l.strip().split()
            if not c[0]=="c:": continue
            if int(c[2])==args.mychunk:
                ccn = int(c[1])
                contigs = eval(" ".join(c[4:]))
                components.append(contigs)
                for cc in contigs:
                    component[cc]=ccn

    if args.progress: print("#Done reading chunkfile")


    promisc = {}
    for n in G.nodes():
        if args.debug: print("#r:",old_div(float(G.degree(n)),ll[n]))
        if (old_div(float(G.degree(n)),ll[n]))>args.promisc: promisc[n]=True

    tag_tallies={}
    bad_tag_tallies={}
    strx={"+":0, "-":1}
    strings = []
    ccn=1
    gn=1

    for c in components:
#    for c in nx.connected_components(SG):
        print(c)
        ccn = component[c[0]]
                
        if args.progress: print("#ccn:",ccn, len(c), sum([ll[cc] for cc in c]))
        gg = nx.subgraph(G,c)
        if args.progress: print("#subgraph with nedges=", len(gg.edges()))
        t= nx.minimum_spanning_tree(gg)
        if args.progress: print("#tree found")
        trim_level=shave_round(t)
        if args.progress: print("#shave round done")
        
        x={}
        for n in t.nodes():
            m = len([ nn for nn in t.neighbors(n) if trim_level[nn]>2 ])
            x[n]=m

        leafDist = distance_to_nearest_leaf(t,trim_level,2,x)
        yDist    = distance_to_nearest_branch(t,trim_level,2,x)

        st = is_shaved_tail(t,trim_level,x,2)

#        ann= annotate_edges(t,G,c)
        bad_edge=False

        for a,b in t.edges():
            aa=tuple()
            bb=tuple()
            if a in besthit: aa = tuple(besthit[a][1:5])
            if b in besthit: bb = tuple(besthit[b][1:5])
            qd = qdist(aa,bb)
            if qd>args.maxLength and (st.get(a,False) and st.get(b,False)) and (t.degree(a)==2 and t.degree(b)==2) and not ( (old_div(float(G.degree(a)), ll[a]))>args.promisc or (old_div(float(G.degree(b)),ll[b]))>args.promisc ) and (trim_level[a]>2) and (trim_level[b]>2) : bad_edge = True

            if (st.get(a,False)==False) and (st.get(b,False)==False) and (not (x[a]>2 and x[b]>2)) and trim_level[a]>2 and trim_level[b]>2: 
                add_tag(a,b,"bigH")

            if args.blacklist and (a,b) in blacklist or (b,a) in blacklist:
                add_tag(a,b,"blacklist")

            if (st.get(a)==True) and (st.get(b)==True): add_tag(a,b,"tail")
            if (promisc.get(a) or promisc.get(b)): add_tag(a,b,"promisc")
 #           if (ann.get((a,b))=='local_bridge'): add_tag(a,b,"bridge")

            if x[a]>2 and x[b]>2 and trim_level[a]>2 and trim_level[b]>2: 
                add_tag(a,b,"H")
                if args.logH: print("#h:",a,b)
            elif (x[a]>2 or x[b]>2) and trim_level[a]>2 and trim_level[b]>2: 
                add_tag(a,b,"Y")
            elif trim_level[a]>2 and trim_level[b]>2 and (yDist.get(a,10) <=2 or yDist.get(b,10)<=2) : 
                add_tag(a,b,"nearY")
#             if trim_level[a]>2 and trim_level[b]>2 and (yDist.get(a,10) <=2 or yDist.get(b,10)<=2) : add_tag(a,b,"nearY")

            if trim_level[a]<=2 or trim_level[b]<=2: 
                add_tag(a,b,"hair")
            elif (trim_level[a]<=3 or trim_level[b]<=3) : 
                add_tag(a,b,"longHair")

            if not args.silent: print("\t".join(map(str,[a,b,ll[a],ll[b],int(-G[a][b]['weight']),G.degree(a),G.degree(b),t.degree(a),t.degree(b),x[a],x[b],trim_level[a],trim_level[b],aa,bb, st.get(a,False), st.get(b,False), leafDist.get(a,-1), leafDist.get(b,-1), yDist.get(a,-1), yDist.get(b,-1) , get_tags(a,b) ])))

            if qd >= args.maxLength : 
                bad_tag_tallies[ get_tags(a,b) ] = bad_tag_tallies.get( get_tags(a,b), 0)+1
#                add_tag(a,b,"good")
#            else:
#                add_tag(a,b,"bad")

            tag_tallies[ get_tags(a,b) ] = tag_tallies.get( get_tags(a,b), 0)+1

            if args.logTags: print("#tags:",a,b,get_tags(a,b))

#        for n in x.keys():
#            if x[n]>2:
#                print "#y:",n,x[n],[ nn for nn in t.neighbors(n) if trim_level[nn]>2 ]

#        d=nx.all_pairs_shortest_path_length(t,cutoff=10)
#        for kk in d.keys():
#            for kkk in d[kk].keys():
#                print kk,kkk,d[kk][kkk]


               # def printdot(g,c,n  ,ll,bh     ,annot,trim_level={},post_trim_degree={},tag="bad",yDist={},leafDist={},edgeTags=edge_tags):
        if bad_edge or args.savetreedots: printdot(t,gg,c,ccn,ll,besthit,{}  , trim_level=trim_level,post_trim_degree=x,tag=args.dotLabel,yDist=yDist,leafDist=leafDist,edgeTags=edge_tags)

        to_remove=[]
        for a,b in t.edges():
            ta=get_tags(a,b)
            if ("Y" in ta) or ("H" in ta) or ("hair" in ta):
                to_remove.append((a,b))
                print(a,b,ta)
        print(to_remove)
        t.remove_edges_from(to_remove)
        seen={}
        backbones=[]
        
        for n in t.nodes():
            if t.degree(n)==1 and not n in seen:
                bb=list(nx.dfs_preorder_nodes(t,n))
                for b in bb: seen[b]=1
                backbones.append( bb )

        singletons=[ n for n in t.nodes() if not n in seen]
#        print "lin:",backbones,singletons

        group={}
        for s in singletons:
            group[s]=0

        for b0 in backbones:
            b,strands=loo.optimize_unoriented_ordering(b0,ll,links,w=3,scaffoldid="{}.{}".format(args.mychunk,gn))
            gn+=1
            print(b,strands)
#            final,score_gap=orient.orient_string(b,1,links,{})
#            print b
#            print final
            for i in range(1,len(b)):
                suff1 = ".3" if strands[i-1]=="+" else ".5" 
                suff2 = ".5" if strands[i  ]=="+" else ".3" 
            
                print("#join:","\t".join(map(str,[b[i-1]+suff1,b[i]+suff2,strands[i-1],strands[i],1])))

#        ccn+=1
        if args.head and ccn>args.head:
            break




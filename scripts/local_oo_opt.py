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
import argparse
import random
import gzip
import glob
import local_scramble_dp as lsd
import chicago_edge_scores as ces
import math
import edgelist2oo as edgelist2oo
import hashlib
import struct
import colorsys
import inspect

default_windowsize=3
default_gapsize=100
global_windowsize=default_windowsize
bl=lsd.LSbacklinks(default_windowsize)
bf=lsd.LSbitlinks(default_windowsize)
tf2pm={True:"-",False:"+"}

def setup_for_windowsize(w):
    global bl
    global bf
    bl = lsd.LSbacklinks(w)
    bf = lsd.LSbitlinks(w)

def optimize_unoriented_ordering(s,ll,links,w=3,scaffoldid=1):
    global global_windowsize
    if not w==global_windowsize:
        global_windowsize=w
        setup_for_windowsize(w)
    r=optimize_oo(s,scaffoldid,ll,w,{},{},dict( [  (x,scaffoldid) for x in s  ]  ), dict( [  (x,"+") for x in s  ]  ) ,links)
    return r

def set_fit_params(fit_params):
    ces.set_exp_insert_size_dist_fit_params(fit_params)

def optimize_oo(sc,scn,ll,w,coords,besthit,scaffold,contig_strand,links,debug=False,viz=False,no_mean_field=True):
        gap_length=default_gapsize
        
        global bl
        global bf
        L=len(sc)

        if debug:
            print("sc:",sc)
            print("scn:",scn)
            print("coords:",coords)
            print("besthit:",besthit)
            print("scaffold:",scaffold)
            print("contig_strand:", contig_strand)
        
        slen=(len(sc)-1)*gap_length
        for s in sc:
            slen += ll[s]

        if False: #L < w : #continue
            plines=[]
            for sca in sc:            
                print("x:",scaffold[sc[0]],0,sca, (0,0), contig_strand[sca]=="-",sc)
#                print "x:",scaffold[sc[0]],0,sca, (0,0), contig_strand[sca]=="-"
                plines.append( [scaffold[sc[0]],sca+".5",int(coords[sca+".5"]),-1,-1,-1,ll[sca],  besthit.get(sca,[-1]*5) ])
                plines.append( [scaffold[sc[0]],sca+".3",int(coords[sca+".3"]),-1,-1,-1,ll[sca],  besthit.get(sca,[-1]*5) ])
            plines.sort(key=lambda x: (x[0],x[2]))
            scl=max( [ p[2] for p in plines] )
            for p in plines:
                p[5]=scl
                print("\t".join(map(str,["p:"]+p)))

            return #continue
        color={}
 

        if debug: print("#",L)
        for s1 in sc:
            color[s1] = colorsys.hls_to_rgb( random.random(), 0.5, 0.5*random.random() + 0.5)
#        for i in range(L-5):
#            print "\t".join(map(str,[sc[i]]+ [ len(links.get((sc[i],sc[j]),[])) for j in range(i+1,i+5) ] ))

        if viz: of=open(viz,"w")
        cs={}

        #flag boundary-violating permutations at the left boundary:
        boundary_crossers={}
        for i in range(-w,0):
            for p in range(bl.n):
                for k in range(0,i+w):
                    if bl.perms[p][k-i] < -i:
                        boundary_crossers[i,p]=True
                        if debug: print("bc:",i,p,k,bl.perms[p][k-i],bl.perms[p])

        #flag boundary-violating permutations at the right boundary:
        for i in range(L-w,L):
            for p in range(bl.n):
                for k in range(i,L):
                    if i+bl.perms[p][k-i] >= L:
                        boundary_crossers[i,p]=True
                        if debug: print("bc:",i,p,k,bl.perms[p][k-i],bl.perms[p],L-i)
                        
        #work out the coordinate system for each scaffold at each DP matrix cell:
        cs={}
        rs=0
        line=0
#        print "cs: --- "
#        print "cs:"
        for i in range(-w+1,L):
            for p in range( bl.n):
                if (i,p) in boundary_crossers: continue
                for f in range( bf.n):
                    line -= 5 
                    lrs=0
                    for k in range(w):
                        flipped=1&(f>>(w-k-1))==1
                        #print i,p,f,k,bl.perms[p],sc[i+bl.perms[p][k]],rs,lrs,flipped
                        if i+bl.perms[p][k]<0 or i + bl.perms[p][k]>=L: continue
                        s1 = sc[i+bl.perms[p][k]]
                        if flipped:
                            cs[i,p,f,k] = (rs+lrs+ll.get(s1,0) ,-1)
                            if viz: of.write("draw(({0},{1})--({2},{1}),rgb({3},{4},{5}),Arrow);\n".format(old_div((rs+lrs+ll.get(s1,0) + gap_length),1000.0), line, old_div((rs+lrs+gap_length),1000.0), color[s1][0],color[s1][1],color[s1][2]))
                        else:
                            cs[i,p,f,k] = (rs+lrs              , 1)
                            if viz: of.write("draw(({0},{1})--({2},{1}),rgb({3},{4},{5}),Arrow);\n".format(old_div((rs+lrs+gap_length),1000.0),line,old_div((rs+lrs+gap_length+ll.get(s1,0)),1000.0),color[s1][0],color[s1][1],color[s1][2]))
#                        if True: print "cs:",i,p,f,k,bl.perms[p],sc[i+bl.perms[p][k]],rs,lrs,flipped,cs[i,p,f,k]
                        lrs+=ll.get( s1,0)+gap_length
            if i>=0:  rs+=ll.get(sc[i],0)+gap_length
        if viz: of.close()

        #sc0

        #initialize the reference scores
        n0   = {}
        sc0  = {}
        cs0  = {}
        rs   = 0
        line = 0

        sc_cache={}

        def get_scpair_score(s1,s2,xdelta,strand1,strand2):
            if (s1,s2,xdelta,strand1,strand2) not in sc_cache:
                #x=0.0
                #for a,b in links.get((s1,s2),()):
                    #x,y =  cs[i,p,f,k][0]+cs[i,p,f,k][1]*a, cs[j,q,g,l][0]+cs[j,q,g,l][1]*b
                #    x += get_score(abs( xdelta + strand1*a - strand2*b ))
                #sc_cache[ s1,s2,xdelta,strand1,strand2 ] = x
                #    def ll(self,l1,l2,o1,o2,links,g,p0=-1):                

#    if (o1,o2) == (0,0):     #         -------1----->    ----2------>
#    if (o1,o2) == (0,1):     #         -------1----->   <----2------
#    if (o1,o2) == (1,0):     #        <-------1-----     ----2------>
#    if (o1,o2) == (1,1):     #        <-------1-----    <----2------

                gaplen = abs(xdelta)
                o1=0
                o2=0
                if xdelta>0: # s1 on the right
                    if strand1==-1: gaplen -= ll[s1]
                    if strand2== 1: gaplen -= ll[s2]

                    if strand1==1: o1=1
                    if strand2==1: o2=1

                else:        # s2 on the right
                    if strand1== 1: gaplen -= ll[s1]
                    if strand2==-1: gaplen -= ll[s2]

                    if strand1==-1: o1=1
                    if strand2==-1: o2=1
                #print(s1,s2)
                try:
                    #ssc = ces.model.ll(ll[s1],ll[s2],o1,o2,links.get((s1,s2),[]),gaplen)
                    ssc = ces.model.score(ll[s1],ll[s2],o1,o2,links.get((s1,s2),[]),gaplen)

                except ValueError as e:
                    print(s1,s2,ll[s1],ll[s2])
                    raise e
                sc_cache[ s1,s2,xdelta,strand1,strand2 ] = ssc
                #print "xy:",s1,s2,strand1,strand2,xdelta,sc_cache[ s1,s2,xdelta,strand1,strand2 ], ssc
            return sc_cache[ s1,s2,xdelta,strand1,strand2 ]

        if not no_mean_field:
            for i in range(-w,L):
                j=i+w
                for p in range(bl.n):                
                    if (i,p) in boundary_crossers: continue
                    for f in range(bf.n):

                        for k in range(w):
                            if i + bl.perms[p][k] < 0 or i + bl.perms[p][k] >=L : continue
                            s1 = sc[ i + bl.perms[p][k] ]

                            for q in range(bl.n):                
                                if (j,q) in boundary_crossers: continue

                                for g in range(bf.n):
                                    for l in range(w):
                                        if j+l -(i +k)>w:
                                            if j + bl.perms[q][l] >=L: continue
                                            s2 = sc[ j + bl.perms[q][l] ]
                                            if debug: print("sc0",i+k,j+l,i + bl.perms[p][k], j + bl.perms[q][l],s1,s2)
                                            n0[s1,s2]= n0.get((s1,s2),0)+1
                                            sc0[s1,s2] = sc0.get((s1,s2),[]) + [get_scpair_score(s1,s2,(cs[i,p,f,k][0]-cs[j,q,g,l][0]),cs[i,p,f,k][1],cs[j,q,g,l][1])]
    #                                        sc0[s1,s2] = sc0.get((s1,s2),0.0) + get_scpair_score(s1,s2,(cs[i,p,f,k][0]-cs[j,q,g,l][0]),cs[i,p,f,k][1],cs[j,q,g,l][1])

                                            #for a,b in links.get((s1,s2),()):
                                            #    x,y =  cs[i,p,f,k][0]+cs[i,p,f,k][1]*a, cs[j,q,g,l][0]+cs[j,q,g,l][1]*b
                                            #    sc0[s1,s2] = sc0.get((s1,s2),0.0) + get_score(abs(y-x))

        setx={}
        for s1,s2 in list(sc0.keys()):
            if debug: print("#d:",s1,s2,n0[s1,s2],sc0[s1,s2],sc0.get((s2,s1),0.0),sc0[s1,s2],n0[s1,s2])
            if no_mean_field:
                sc0[s1,s2]=0.0
                sc0[s2,s1]=0.0
            else:
                if (s1,s2) in setx or (s2,s1) in setx:  
                    print("double hit!")
                    exit(0)
                setx[s1,s2]=1
                setx[s2,s1]=1
                sc0[s1,s2] = max(sc0[s1,s2])
                sc0[s2,s1] = sc0[s1,s2]

        S={}
        bt={}
        minus_infty=-9.9e99

        #set the boundary conditions:
        for i in range(-w,0):
            for p in range(bl.n):
                
                for f in range(bf.n):
                    score = 0.0
                    if (i,p) in boundary_crossers: 
                        #S[0,p,f] = minus_infty
                        continue

                    for k in range(w):
                        if i + bl.perms[p][k]< 0 or i + bl.perms[p][k]>=L: continue
                        s1 = sc[ i + bl.perms[p][k] ]

                        for j in range(k+1,w):
                            if i + bl.perms[p][j]< 0 or i + bl.perms[p][j] >=L: continue
                            s2 = sc[ i + bl.perms[p][j] ]
                            lk = links.get(( s1,s2),())

                            score_delta=get_scpair_score(s1,s2,(cs[i,p,f,k][0]-cs[i,p,f,j][0]),cs[i,p,f,k][1],cs[i,p,f,j][1])
                            #for a,b in lk:
                            #x,y =  cs[i,p,f,k][0]+cs[i,p,f,k][1]*a, cs[i,p,f,j][0]+cs[i,p,f,j][1]*b
                                #x,y = cs[0,p,f,i][0]+cs[0,p,f,i][1]*a, cs[0,p,f,j][0]+cs[0,p,f,j][1]*b
                                #if not x<y:
                                #    print i,p,f,j,q,g,s1,s2,a,b,x,y,y-x #,math.log(ces.insert_size_dist(y-x))
                                #score_delta += get_score(abs(y-x)) # math.log(ces.insert_size_dist(y-x))
                            score+=score_delta - sc0.get((s1,s2),0.0)
                            if debug: print("score delta", s1,s2,i,p,f,score_delta)
                    S[i,p,f] = score
                
        #fill in the matrix:
        for i in range(0,max(L-w+1,1)):
            if debug: print("#i:",i)
            for p in range( bl.n):
                if debug: print("#p:",p)
                for f in range( bf.n):
                    if debug: print("#f:",f)
                    best_score=-1.0e9
                    this_backtrace_step=()
                    #                    any_back_links=False

                    for d,q in bl.backlinks[p]:
                        j=i-d
                        if debug: print(j,q)
                        if (j,q) in boundary_crossers: continue

                        for g in bf.bitlinks[f,d]:
                            #                            score = S.get((j,q,g),0.0)
                            score = 0.0 #S[j,q,g]
                            #any_back_links=True
                            
                            for dl in range(w-d,w):                                
                                l = i+dl
                                if l<0 or l>=L: continue
                                if i + bl.perms[p][dl]  >= L :continue
                                s2 = sc[ i + bl.perms[p][dl] ]

                                #pairs of contigs where the first is already scored by the back-linked state
                                for dk in range(w):
                                    k = j+dk
                                    if k<0 or k >=L: continue
                                    s1 = sc[ j + bl.perms[q][dk] ]
                                    if (l-k)<=w:
                                        #lk = links.get(( s1,s2),()) 
                                        #print k,l,s1,s2,lk
                                        score_delta=0.0
                                        score_delta=get_scpair_score(s1,s2,(cs[j,q,g,dk][0]-cs[i,p,f,dl][0]),cs[j,q,g,dk][1],cs[i,p,f,dl][1])

                                        #incrs=[]
                                        #for a,b in lk:
                                        #    x,y = cs[j,q,g,dk][0]+cs[j,q,g,dk][1]*a, cs[i,p,f,dl][0]+cs[i,p,f,dl][1]*b
                                        #    if not x<y:
                                        #        pass
                                        #        if args.debug: print i,p,f,j,q,g,s1,s2,a,b,x,y,y-x #,math.log(ces.insert_size_dist(y-x))
                                            #incrs.append(get_score(y-x))
                                        #    score_delta += get_score(y-x) # math.log(ces.insert_size_dist(y-x))
                                        #print "score delta",s1,s2,score_delta, incrs,"#"
                                        score+=score_delta - sc0.get((s1,s2),0.0)
                                        if debug : print("\t".join(map(str,["score delta", score, score_delta, sc0.get((s1,s2)),score_delta - sc0.get((s1,s2),0.0), s1,s2,dk,dl,j,i,i + bl.perms[p][dl],j + bl.perms[q][dk], len(links.get(( s1,s2),())),i,p,f,j,q,g,k,l,l-k])))

                                    else:
                                        #print "skip",l,k,l-k
                                        pass

                                #pairs of contigs where both are being newly scored here
                                for dk in range(w-d,dl):
                                    k = i+dk
                                    if i + bl.perms[p][dk] >=L: continue
                                    s1 = sc[ i + bl.perms[p][dk] ]
                                    if (l-k)<=w:
                                        #lk = links.get(( s1,s2),()) 
                                        #print k,l,s1,s2,lk
                                        #incrs=[]
                                        score_delta=get_scpair_score(s1,s2,(cs[i,p,f,dk][0]-cs[i,p,f,dl][0]),cs[i,p,f,dk][1],cs[i,p,f,dl][1])

                                        #for a,b in lk:
                                        #    x,y = cs[i,p,f,dk][0]+cs[i,p,f,dk][1]*a, cs[i,p,f,dl][0]+cs[i,p,f,dl][1]*b
                                        #    if not x<y:
                                        #        pass
                                        #        if args.debug: print i,p,f,j,q,g,s1,s2,a,b,x,y,y-x #,math.log(ces.insert_size_dist(y-x))
                                                #score += get_score(y-x) # math.log(ces.insert_size_dist(y-x))
                                        #    score_delta += get_score(y-x) # math.log(ces.insert_size_dist(y-x))
                                            #incrs.append( get_score(y-x))
                                        #print "score delta",s1,s2,score_delta, incrs
                                        #score+=score_delta - sc0[s1,s2]
                                        score+=score_delta - sc0.get((s1,s2),0.0)
                                        if debug : print("\t".join(map(str,["score deltaX", score, score_delta, sc0.get((s1,s2)),score_delta - sc0.get((s1,s2),0.0), s1,s2, dk,dl,j,i,i + bl.perms[p][dl],i + bl.perms[p][dk], len(links.get(( s1,s2),())),i,p,f,j,q,g,k,l,l-k])))
                                    else:
                                        #print "skip",l,k,l-k
                                        pass
                            tag=""
                            if score + S[j,q,g]> best_score:
                                best_score = score + S[j,q,g]
                                this_backtrace_step = (j,q,g)
                                tag="XXX"
                            if debug: print("score",i,p,f,j,q,g,S[j,q,g],score , S[j,q,g]+score, bl.perms[p],best_score,tag)

                    S[i,p,f] = best_score
                    bt[i,p,f]=this_backtrace_step
        i=max(0,L-w) # 0 only for scaffolds with len < w
        overall_best_score = minus_infty
        state=()
        for p in range(bl.n):
            if (i,p) in boundary_crossers: continue
            for f in range( bf.n):
                if overall_best_score < S[i,p,f]:
                    overall_best_score = S[i,p,f]
                    state=(i,p,f)

        i,p,f=state
        if debug: print("overall best:",state, overall_best_score, [ sc[i+bl.perms[p][dk]] for dk in range(w)])
        last_i=L

        plines=[]
        my_ordering=[]
        my_orientations=[]
        while state[0]>=-w+1:
            i,p,f = state
            if debug: print("bt:",L,state,state[0],S[state])#-bt[state][0]
            for j in range(last_i-1,max(-1,i-1),-1):
                if i+bl.perms[ p ][j-i] >=L : continue
                flipped=1&(f>>(w-(j-i)-1))==1
                contig = sc[ i+bl.perms[ p ][j-i] ]
                #print i,j,contig
                my_ordering = [ contig] + my_ordering
                my_orientations = [ tf2pm[flipped] ] + my_orientations
                print("x:",scaffold[sc[0]],j,sc[ i+bl.perms[ p ][j-i] ], cs[i,p,f,j-i], flipped) 
                #                if flipped:
                    #                    p3 = 
                x,chr,y = coords.get((contig,"5"),(-1,-1,-1))
                x = cs[i,p,f,j-i][0]
                bh = besthit.get(contig,[-1]*5)
                chr=bh[1]
                if bh[2]=="+":
                    y=int(bh[3])
                else:
                    y=int(bh[4])
#                print "p:",scaffold[contig],contig+".5",x,chr,y,slen[scaffold[contig]],ll[contig],bh
#                print "p:",scaffold[contig],contig+".5",x,chr,y,slen[scaffold[contig]],ll[contig],bh
                plines.append( [scaffold[contig],contig+".5",x,chr,y,slen,ll[contig],bh] )
                x,chr,y = coords.get((contig,"3"),(-1,-1,-1))
                chr=bh[1]
                x = cs[i,p,f,j-i][0] +cs[i,p,f,j-i][1]*ll[contig] 
                if bh[2]=="-":
                    y=int(bh[3])
                else:
                    y=int(bh[4])
#                y=int(bh[4])
#                print "p:",scaffold[contig],contig+".3",x,chr,y,slen[scaffold[contig]],ll[contig],bh
    #                print "p:",scaffold[contig],contig+".3",x,chr,y,slen[scaffold[contig]],ll[contig],bh
                plines.append( [scaffold[contig],contig+".3",x,chr,y,slen,ll[contig],bh]) 
            state = bt.get(state,(-w,0,0))
            last_i=i
        plines.sort(key=lambda x: (x[0],x[2]))
        for p in plines:
            print("\t".join(map(str,["p:"]+p)))
        return( my_ordering,my_orientations)
#        print my_orientations


if __name__=="__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i','--input',default="-")
    parser.add_argument('-X','--mean_field',action="store_true",default=False)
    parser.add_argument('-a','--mychunk',default=0,type=int)
    parser.add_argument('-N','--nchunks',default=1,type=int)
    parser.add_argument('-S','--scaffold',default=[],action="append",help="List of scaffolds to run on")
    parser.add_argument('-l','--links')
    parser.add_argument('-b','--besthits')
    parser.add_argument('-d','--debug',default=False,action="store_true")
    parser.add_argument('-V','--viz',default=False)
    parser.add_argument('-E','--edgelist',action="store_true",default=False)
    parser.add_argument('-w','--window',type=int,default=3)
    parser.add_argument('--seed',required=False,type=int,default=1, help="Seed for random number generation, use -1 for no seed")
    parser.add_argument('-M','--set_insert_size_dist_fit_params') #,default="3.85301461797326,1.42596694138494,1.38674994280385e-05,10940.8191219759,49855.7525034142,0.3,420110993")
    args = parser.parse_args()
    if args.seed != -1 :
        random.seed(args.seed)



    print("#"+str(args))
    if False:
        if args.set_insert_size_dist_fit_params:
            s=args.set_insert_size_dist_fit_params
            a,b,c,d,f,pn,N = list(map(float,list(s.split(','))))
            print("#",a,b,c,d,f)
            ces.set_insert_size_dist_fit_params(a,b,c,d,f,pn,N)
            sys.stdout.flush()

    else:

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
            ces.set_exp_insert_size_dist_fit_params(fit_params)
 

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


    #    parser.add_argument('-s','--scaffold',default=50)
    #parser.add_argument('-o','--outfile')
    #parser.add_argument('-c','--chr',default="3")
    #parser.add_argument('-m','--min',default=0.0   ,type=float)
    #parser.add_argument('-a','--max',default=10.0e6,type=float)


    ll={}
    contig_coord={}
    coords={}
    coord_map={}
    scaffold={}
    contig_strand={}
    bh={}
    if args.input=="-":
        f=sys.stdin
    else:
        f=open(args.input)

    scaffolds=[]
    last_sc=-1
    slen={}
    contig2scn={}
    scaffold_hashes={}
    if args.edgelist:
    #    return(contigs,strand,ll,scaffolds,scaffold,coords)

        scaffolds,contig_strand,ll,scaffold_lists,contig2scn,coords = edgelist2oo.edges2oo(f)
        print("#")
#        for c,s in contig2scn.items():
#            print "contig2scn:",c,s
        scn=1
        for sc in scaffolds:
            sll=0
            for c in sc:
#                print c,ll.get(c)
                sll+=ll.get(c,0)
#                scaffold[c]=scn
                scaffold[c]=contig2scn.get(c,scn)
#                if scaffold[c]==-1: print "missing scaffold info for",c
            slen[scaffold[sc[0]]]=sll
#            slen[scn]=sll

            scn+=1
    else:
        while True:
            l=f.readline()
            if not l: break
            if l[:2]=="p:":
                c=l.strip().split()
                #print c
                v=eval(" ".join(c[8:]))
                strand="+"
                if v:
                    strand = v[2]
                #sc,contig_end,x,chr,y,sslen = int(c[1]),c[2],float(c[3]),c[4],float(c[5]),float(c[6])#,int(c[7])
                sc,contig_end,x,chr,y,sslen,llen = c[1],c[2],float(c[3]),c[4],float(c[5]),float(c[6]),int(c[7])
                if not sc == last_sc:
                    if len(scaffolds)>0 and len(scaffolds[-1])==0: sys.stderr.write("wtf? empty list for scaffold {} at {}\n".format(last_sc,inspect.getframeinfo(inspect.currentframe())[:3] ))
                    scaffolds.append([])
                    #print sc,last_sc,scaffolds

                contig,end = contig_end[:-2],contig_end[-1:]
                ll[contig]=llen

                if end=="3":
                    scaffolds[-1].append(contig)

                contig_strand[contig]=strand

                coords[contig,end]=(x,chr,y)
                coords[contig_end]=x #(x,chr,y)
                slen[sc]=sslen
                scaffold[contig]=sc

                if not sc in scaffold_hashes:
                    h=struct.unpack("<L", hashlib.md5(sc.encode("utf-8")).digest()[:4])[0]
                    scaffold_hashes[sc]=h%args.nchunks 
                    #print("hash:",sc,scaffold_hashes[sc])

                bh[contig]=v
                last_sc=sc

        for i in range(len(scaffolds)):
            scaffolds[i].sort( key= lambda x: coords[x,"3"] )

    links={}
    if args.links:
        for fn in glob.glob(args.links):
            print("#",fn)
            if fn[-3:]==".gz":
                f = gzip.open(fn)
            else:
                f = open(fn)

            while True:
                l=f.readline()
                if not l: break
                if l[0]=="#": continue
                c=l.strip().split()
                if not ( c[0] in scaffold and c[1] in scaffold and scaffold[c[0]]==scaffold[c[1]] ):
                    continue
                v=eval(" ".join(c[5:]))
                s1,s2,l1,l2,nl = c[0],c[1],int(c[2]),int(c[3]),int(c[4])
                links[s1,s2]=v
                #print l.strip(),"#",scaffold[c[0]]

    w=args.window
    if not w==default_windowsize:
        setup_for_windowsize(w)
    gap_length=default_gapsize

    for sc in scaffolds:
        #print("z",sc[0],scaffold[sc[0]],scaffold_hashes[scaffold[sc[0]]],args.mychunk,sc[:3])
        if len(sc)==0: continue
        if ((not args.edgelist) and (not scaffold_hashes[scaffold[sc[0]]]==args.mychunk)) or (args.scaffold and not scaffold[sc[0]] in args.scaffold): continue

#        for sca in sc: print sca
        scn = contig2scn.get(sc[0],False)
        L=len(sc)
        if scn: print("cc:",scn,L,scaffold_lists.get(scn,[]))

        r=optimize_oo(sc,scn,ll,w,coords,besthit,scaffold,contig_strand,links,debug=args.debug,viz=args.viz,no_mean_field=not args.mean_field)
        print(r)

#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
import sys
import argparse
from mapperlite import MapperLite
import struct
import hashlib
import pysam
import chicago_edge_scores as ces
#import BamTags
from bamtags import BamTags
from chicago_edge_links2 import read_mask_ranges, mask_test
from time import gmtime, strftime
debug=False

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

def segments_intersect(x,y):
    if not pairs_overlap(x,y):
        raise Exception 

    a=min(x[0],x[1])
    b=max(x[0],x[1])
    c=min(y[0],y[1])
    d=max(y[0],y[1])
        
    if a<=c and c<=b: 
        p,q = (c,min(b,d))
        if p>q: 
            print("wtf1",x,y,(p,q))
    if a<=d and d<=b: 
        p,q = (max(a,c),d)
        if p>q: 
            print("wtf2",x,y,(p,q))
    if c<=a and a<=d: 
        p,q = (a,min(b,d))
#wtf3 (2608, 2741) (1500, 3000) (2608, 1500)
        if p>q: 
            print("wtf3",x,y,(p,q))
    if c<=b and b<=d: 
        p,q = (max(a,c),b)
        if p>q: 
            print("wtf4",x,y,(p,q))

#wtf (1694, 3362) (1500, 3000) (1694, 1500)

    if p>q: 
        print("wtf",x,y,(p,q))
        raise Exception
    return (p,q)

tmsd_debug=False

def tile_mask_score_delta(i,j,binwidth,masked_segments,model,mask_iter_i,mask_iter_j,debug):
    if i==j: return(0)
    gap=(j-i)*binwidth
    left_tile_masked_segments  = []
    right_tile_masked_segments = []
    mii0,mij0=mask_iter_i,mask_iter_j
    while mask_iter_i<len(masked_segments) and masked_segments[mask_iter_i][0]<(i+1)*binwidth:
        if tmsd_debug: print("tmsd:\tI",mask_iter_i,masked_segments[mask_iter_i],(i*binwidth,(i+1)*binwidth),file=sys.stderr,sep="\t")
        if pairs_overlap(            masked_segments[mask_iter_i],(i*binwidth,(i+1)*binwidth)) : 
            a,b = segments_intersect(masked_segments[mask_iter_i],(i*binwidth,(i+1)*binwidth))
            if b-a < 0:
                print("tmsd:\tI",a,b,mask_iter_i,masked_segments[mask_iter_i],(i*binwidth,(i+1)*binwidth),file=sys.stderr,sep="\t")
                raise Exception
            left_tile_masked_segments.append(  (a,b) )
        mask_iter_i+=1

    while mask_iter_j<len(masked_segments) and masked_segments[mask_iter_j][0]<(j+1)*binwidth:
        if tmsd_debug: print("tmsd:\tJ",mask_iter_j,masked_segments[mask_iter_j],(j*binwidth,(j+1)*binwidth),file=sys.stderr,sep="\t")
        if pairs_overlap(            masked_segments[mask_iter_j],(j*binwidth,(j+1)*binwidth)) : 
            a,b = segments_intersect(masked_segments[mask_iter_j],(j*binwidth,(j+1)*binwidth))
            if b-a < 0:
                print("tmsd:\tJ",a,b,mask_iter_j,masked_segments[mask_iter_j],(j*binwidth,(j+1)*binwidth),file=sys.stderr,sep="\t")
                raise Exception
            right_tile_masked_segments.append( (a,b) )
        mask_iter_j+=1

    score_delta = 0.0

    for a,b in right_tile_masked_segments:
        score_delta += model.n_bar( binwidth, b-a, gap+a-j*binwidth ) - model.n_bar0( binwidth*(b-a) )

    for a,b in left_tile_masked_segments:
        score_delta += model.n_bar( binwidth, b-a, gap+(i+1)*binwidth-b ) - model.n_bar0( binwidth*(b-a) )

    for a,b in left_tile_masked_segments:
        for c,d in right_tile_masked_segments:
            if tmsd_debug: print("mask pair:",a,b,c,d,b-a,d-c,c-b,file=sys.stderr,sep="\t")
            score_delta -= model.n_bar( b-a, d-c, c-b ) - model.n_bar0( (b-a)*(d-c) )

    if tmsd_debug:
        print("tmsd:",(i*binwidth,(i+1)*binwidth),(j*binwidth,(j+1)*binwidth),score_delta,left_tile_masked_segments,right_tile_masked_segments,(mii0,mij0),masked_segments,file=sys.stderr,sep="\t")
    
    return(score_delta)

def pair2bin2D(pair,binwidth):
     a,b=pair
     return( int(a/binwidth),int(b/binwidth) )

def make_scaffold_mask(scaffold,mapper,mask):
    s_mask = {}
    segments={}
    slen = mapper.scaffold_length[scaffold]
    for seg in mapper.scaffold_contigs[scaffold]:
#        print("x:",seg)
        ocontig=mapper.contig_ocontig[seg]
        for m in mask.get(ocontig,[]):
            oscaffold ,z1a,z2a,strand,c1 = mapper.mapCoord(ocontig,m[0],m[1])
            if not scaffold==oscaffold: continue
            if not z1a: continue
            if strand==1:
#                segments.append( (z1a,min(slen,z2a)) )
                segments[  (z1a,min(slen,z2a))  ] =1      
            else:
                segments[ (max(0,z2a),z1a) ] = 1
#                segments.append( (max(0,z2a),z1a) )

#            print(ocontig,seg,"\t\t",m,mapper.mapCoord(ocontig,m[0],m[1]))
#    for a,b in segments:
#        print("\t",a,b,b-a,sep="\t")
    segments=list(segments.keys())
    segments.sort()
    return(segments)

def chicago_pairs(sca,mapper,bamlist,minq=20,mask={}):
    for seg in mapper.scaffold_contigs[sca]:
          ref="_".join(seg.split("_")[:-1])
          for b in bamlist:
               for aln in b.fetch(until_eof=True,reference=ref):
                    if not aln.is_read1: continue
                    if aln.is_duplicate: continue
                    if aln.mapq < minq : continue
                    if BamTags.mate_mapq(aln) < minq : continue
#                    print("#x",ref,mask.get(ref,[]))
                    if mask_test(ref,aln.pos,mask) or mask_test(b.getrname(aln.rnext),aln.pnext,mask) : continue
                    
                    contig = b.getrname(aln.tid) # snam[aln.tid]
                    
                    ncontig= b.getrname(aln.rnext) if aln.rnext>=0 else -1

                    scaffold ,z1a,z2a,z3a,c1 = mapper.mapCoord( contig, aln.pos,     aln.pos+1 ) 
                    nscaffold,z2p,x2p,z3p,c2 = mapper.mapCoord(ncontig, aln.pnext,aln.pnext+1 ) 
                    if debug: print(("#x",contig,ncontig,aln.pos,aln.pnext,scaffold,nscaffold,sca,z1a,z2p,ref,mapper.ocontig_contigs.get(ref,[])))
                    if scaffold==nscaffold and sca==scaffold:
                         #yield( sc,seg,contig,ncontig,scaffold,z1a,z2a,z3a,nscaffold,z2p,x2p,z3p )
                         yield( sca,z1a,z2p,c1,c2,seg,contig,ncontig,scaffold,z1a,z2a,z3a,nscaffold,z2p,x2p,z3p,aln.query_name )


hist_sample_interval = 1500
hist_bisize = 1.0
scaffold_end_filter=10000
# Fine grained means for every read
def fine_grained_support(edges,nb,scaffold,pairs,model,buff,debug,t1,t2,gpf,minx,maxx,joins,logfile,binwidth,slen,mask=[],masked_segment_pairs_n_minus_n0=[],raw_hist=False,support_curve=False):

#     print("zz:",scaffold,sum([b-a for a,b in mask]),[ (a,b,b-a) for a,b in mask])

     curve=[]
     if gpf: gpf.write("\n\n\n")

#     if  True: print("fine",len(edges),debug)
     if len(edges)==0: return([],[])
     if edges[-1][0]<0 or edges[-1][0] > slen+1: 
         print("unexpected coordinate {} not in range 0-{}".format(edges[-1][0],slen))
         raise Exception 

     edges.sort()
     nb.sort()
     rs=0 # running sum for individual reads
     n=0
     tripped=False # have we reached the threshold where we want to start breaking (by construction support low on edges).
     last=0
     state=0
     stretch_start=0
     low_point=0.0
     ji=0
     gap_scores={}
     gap_coverages={}
#     slen = mapper.scaffold_length[scaffold]
     break_buffer = []
     last_trip=0

     a=0 
     b=0
     c=0
     f=0 # index in masked segment pairs

     mask_correction_rs=0.0
     next_sample = hist_sample_interval /2
     for i in range(len(edges)):
          rs+=edges[i][1]
          n+=nb[i][1] # number of pairs spanning
          x=edges[i][0] # position

          while f<len(masked_segment_pairs_n_minus_n0) and masked_segment_pairs_n_minus_n0[f][0]<x: 
              mask_correction_rs -= masked_segment_pairs_n_minus_n0[f][1]
              f +=1

          #update iterators that keep track of which masked regions are "in range":  a-b are in range and trailing x, b-c are in range and upcoming.
          while a<len(mask) and mask[a][1]<x-maxx: a+=1
          while b<len(mask) and mask[b][0]<x:      b+=1
          while c<len(mask) and mask[c][0]<x+maxx: c+=1
          # -----------------------------------------
          # a               bx                      c
          mask_correction = 0.0
          # rmi = right mask index
          for rmi in range(b,min(c+1,len(mask))):
              ma,mb = mask[rmi]
              if x>ma: continue
              left_limit = max(0,ma-maxx)
              if (x-left_limit)<0: continue
              mask_correction_delta = model.n_bar( mb-ma, x-left_limit, ma-x ) - model.n_bar0( (mb-ma)*(x-left_limit ) )
              mask_correction += mask_correction_delta
              if debug: print(x,ma-x,a,b,c,left_limit,"Right",mask_correction_delta,sep="\t")

          # left mask index
          for lmi in range(a,b):
              ma,mb = mask[lmi]
              if x<mb:
                  if ma<x and x<mb: 
                      right_limit = min(slen,mb+maxx)
                      #left_limit = min(slen,mb-maxx)
                      left_limit = max(0,ma-maxx)
                      mask_correction_delta1 = model.n_bar( mb-x, x-left_limit, 0 )     - model.n_bar0( (mb-x)*(x-left_limit) )
                      try:
                          mask_correction_delta2 = model.n_bar( x-ma, right_limit-x, mb-x ) - model.n_bar0( (x-ma)*(right_limit-x) )
                      except Exception as e:
                          #wtf: scaffold: 7396, x: 1006292, ma: 686903, mb: 1006300, rl: 886903, ll: 486903
                          print("wtf: scaffold: {scaffold}, x: {x}, ma: {ma}, mb: {mb}, rl: {right_limit}, ll: {left_limit}".format(scaffold=scaffold,x=x,ma=ma,right_limit=right_limit,mb=mb,left_limit=left_limit))
                          raise e
                      mask_correction += mask_correction_delta1 + mask_correction_delta2          
                      if debug: print(x,x-mb,a,b,c,left_limit,right_limit,"Spanned",mask_correction_delta1,mask_correction_delta2,sep="\t")
                      
                  continue
              right_limit = min(slen,mb+maxx)
              if right_limit - x < 0: continue                  
              mask_correction_delta = model.n_bar( mb-ma, right_limit-x, x-mb ) - model.n_bar0( (mb-ma)*(right_limit-x ) )
              mask_correction += mask_correction_delta          
              if debug: print(x,x-mb,a,b,c,right_limit,"Left",mask_correction_delta,sep="\t")

          # take our running sum score and see if we should cut here, fine grained not clipping
          try:
               score=model.cutScore(slen,x,n,rs,rangeCutoff=maxx,minSep=minx)+mask_correction+mask_correction_rs
               score+=mask_correction+mask_correction_rs
               if debug: print("finescore:",scaffold,x,a,b,c,len(mask),score,mask_correction ,mask_correction_rs,score+mask_correction+mask_correction_rs,sep="\t")
          except Exception as e:
               print("Exception computing cutScore for scaffold {} at {}. i={}".format(scaffold,x,i),edges[:10])
               #Exception computing cutScore for scaffold 8351 at 12290. i=2 [(11306, -8.997670875606678, 'a', 'a'), (11686, -8.578118407318865, 'a', 'a'), (12290, 8.578118407318865, 'a', 'b'), (12297, 8.997670875606678, 'a', 'b')]
#               print(i,edges[i],nb[i],rs,n,x,scaffold,slen,len(edges))
               raise e

#          print("#xxxxx",x,next_sample,raw_hist)
#          print(support_curve)
          if support_curve: 
              curve.append((x,score))
#              print("##yield")
#              yield (x,score) 
          if (not raw_hist==False) and x>next_sample:
              if min(x,slen-x)>scaffold_end_filter:
                  next_sample += hist_sample_interval
                  hist_bin=int(score/hist_bisize)*hist_bisize
                  raw_hist[ hist_bin ] = raw_hist.get(hist_bin,0)+1
#              print("#hist bin",hist_bin,score,raw_hist[ hist_bin ])
#              sys.stdout.flush()

          # Obsolete? each gap between contigs?
          while ji<len(joins) and x>joins[ji][0]: 
               gap_scores[ji]    =score
               gap_coverages[ji] =n
               ji+=1

          if score>t1:
               tripped=True
               last_trip=x
          if tripped and score<t2 and state==0:
               stretch_start=x
               state=1
               low_point =score
               low_x = x

          # Reached beginning of region for candidate break (state = 0)
          if state==1 and score>t2:
               
               if logfile: 
                    break_buffer.append( (scaffold,stretch_start,x,low_x,low_point,slen)  )
#               print(scaffold,stretch_start,x,slen,low_point,"break")
               state=0
          if state==1:
               if score<low_point:
                    low_point=score
                    low_x=x
          if debug: print("dd:",scaffold,edges[i][0],score,rs,state,x,stretch_start,score,n,edges[i][2])
          if gpf: gpf.write("{}\t{}\t{}\t{}\n".format(edges[i][0],score,rs,n))

          last=edges[i][0]

     segments = []
     for scaffold,stretch_start,x,low_x,low_point,slen in break_buffer:
          if x < last_trip:
               logfile.write("{} {} {} {} {} {} rawLLR\n".format(scaffold,stretch_start,x,low_x,low_point,slen))
               segments.append((scaffold,stretch_start,x,low_x,low_point))

     return segments,curve
#     for ji in range(len(joins)):
#          print("\t".join(map(str,["gapscore:",ji]+list(joins[ji])+[gap_scores.get(ji,0),gap_coverages.get(ji,0)])))

def pairs2support(scaffold,
                  pairs, # internal to scaffold
                  model,
                  slen=0,
                  masked_segments=[],
                  mapper=None,
                  buff=0,
                  debug=False,
                  t1=20.0,
                  t2=5.0,
                  gpf=False, # gnu plot file handle?
                  minx=1000,
                  maxx=1e7,
                  joins=[],
                  logfile=False,
                  binwidth=1000,
                  nreject=2, # how many rows or colums to toss out
                  raw_hist=False, 
                  clipped_hist=False,
                  support_curve=False):


    
#     slen=mapper.scaffold_length[scaffold]
#     if debug: print("#masked:",scaffold,len(masked_segments),masked_segments[:10]) 
#     print("#",scaffold,masked_segments,file=logfile,sep="\t")
     logfile.flush()

     edges=[] # pairs
     nb=[] # buffer to keep track of how many pairs cover
     
     tile_scores={}
     tile_counts={}
     tile_reads={}
     maxtile=0
#     masked_segments=[]
     for p in pairs:
          if len(p)>=16: # old style, mostly ignore the additional info now
              if p[1]<p[2]:
                   a,b,c1,c2,w,z=p[1],p[2],p[3],p[4],p[6],p[7]
              else:
                   a,b,c1,c2,w,z=p[2],p[1],p[4],p[3],p[7],p[6]
    #          masked_segments = p[17]
              tid=p[16]
          else:
#              a,b = p[0],p[1]
              # Order the coordinates
              a=min( p[0],p[1])
              b=max( p[0],p[1])
              c1,c2="a","a"
              tid="x"

          if a>slen or b> slen:
              print("how could a read be at x > slen?",a,b,slen)
              raise Exception 

          # Throw away really far and short innies
          if abs(b-a) > maxx: continue
          if abs(b-a) < minx: continue
          
          # insert size log liklihood
          ll=model.lnF(b-a)

          if debug: print("pt:",scaffold,a,b,b-a,tid,ll-model.lnFinf)

          #For the old-style exact LLR score, size of "buff" is taken off
          # to be conservative about support
          edges.append( tuple([a+buff         ,  ll,c1,"a"]) )
          edges.append( tuple([b-buff +1      , -ll,c2,"b"]) )
          nb.append( tuple([a+buff         , 1]) )
          nb.append( tuple([b-buff    +1   ,-1]) )

          #For the new-style clipped LLR score add to appropriate tile
          tile = pair2bin2D((a,b),binwidth) # tile is a tuple rows, colums
          maxtile = max(maxtile,tile[0],tile[1]) # furthest tiles seen for size?
          tile_scores[tile] = tile_scores.get(tile,0.0)+ll
          tile_counts[tile] = tile_counts.get(tile,0)+1
          # for debuggin? should we remove when not debugging?
          tile_reads[tile] = tile_reads.get(tile,[]) + [tid]

     masked_segment_pairs_n_minus_n0 = []
     for i in range(len(masked_segments)):
         a,b = masked_segments[i]
         for j in range(i+1,len(masked_segments)):
             c,d = masked_segments[j]               #   pair of masked segments:   a---b     c-----d
             if c-b > maxx: continue                #   gap = c-b
             if b-a<0: 
                 print("#wtf?",scaffold,(a,b),(c,d),i,j,file=logfile,sep="\t")
                 logfile.flush()
                 raise Exception 
             if d-c<0: 
                 print("#wtf?",scaffold,(a,b),(c,d),i,j,file=logfile,sep="\t")
                 logfile.flush()
                 raise Exception 
             if c-b<0: 
                 print("#wtf?",scaffold,(a,b),(c,d),i,j,file=logfile,sep="\t")
                 logfile.flush()
                 raise Exception 
             
             n_bar = ces.model.n_bar(b-a,d-c,c-b)
             n_bar0= ces.model.n_bar0((b-a)*(d-c))
             if debug: print("X",i,j,(a,b),(c,d),b-a,d-c,c-b,n_bar,n_bar0,sep="\t")
             masked_segment_pairs_n_minus_n0.append(  (b, (n_bar-n_bar0))  )
             masked_segment_pairs_n_minus_n0.append(  (c,-(n_bar-n_bar0))  )
     masked_segment_pairs_n_minus_n0.sort()
     fine_grain_segments,fine_grain_support_curve=fine_grained_support(edges,nb,scaffold,pairs,model,buff,debug,t1,t2,gpf,minx,maxx,joins,logfile,binwidth,slen,masked_segments,masked_segment_pairs_n_minus_n0=masked_segment_pairs_n_minus_n0,raw_hist=raw_hist,support_curve=support_curve)

     if debug:
         print("w:",scaffold,slen,len(fine_grain_segments),len(edges),sep="\t")
         for a,b in masked_segments:
             print(a,b)
             if a>slen: 
                 print("a>slen",a,slen)
                 raise Exception 
             if b>slen: 
                 print("b>slen",b,slen)
                 raise Exception 

     # clipped LLR score:
     tile_bumps=[]


#     print("#maxtile=",maxtile,strftime("%Y-%m-%d %H:%M:%S"),file=logfile,sep="\t")
#     logfile.flush()
     mask_iter_i=0
     for i in range(maxtile+1):
#          print("#i=",i,strftime("%Y-%m-%d %H:%M:%S"),file=logfile,sep="\t")
#          logfile.flush()

          while mask_iter_i<len(masked_segments) and masked_segments[mask_iter_i][1]<i*binwidth: mask_iter_i+=1
          j=i          
          mask_iter_j=mask_iter_i
          while ((j-i)*binwidth<maxx) and j<=maxtile:
               tile=(i,j)
#               print("#tile=",(i,j),strftime("%Y-%m-%d %H:%M:%S"),file=logfile,sep="\t")
#               logfile.flush()
               tscore = tile_scores.get(tile,0.0)
               score = model.tileScore(binwidth,tile,tile_counts.get(tile,0),tscore,rangeCutoff=maxx,minSep=minx)
#               print("#score=",score,masked_segments,i,j,mask_iter_i,mask_iter_j,strftime("%Y-%m-%d %H:%M:%S"),file=logfile,sep="\t")
#               logfile.flush()

               mask_score_delta = tile_mask_score_delta(i,j,binwidth,masked_segments,model,mask_iter_i,mask_iter_j,debug)
               score+= mask_score_delta



#               print("#mask_score_delta=",mask_score_delta,strftime("%Y-%m-%d %H:%M:%S"),file=logfile,sep="\t")
#               logfile.flush()

               if debug: 
                    print("tile:",scaffold,tile[0],tile[1],tile[0]*binwidth,tile[1]*binwidth,tscore,tile_counts.get(tile,0.0),score,mask_score_delta,sep="\t")
                    for read in tile_reads.get(tile,[]):
                         print("tileread:",tile,read)
               if not i==j:
                    tile_bumps.append( (i*binwidth,i,j, score, 1)  )
                    tile_bumps.append( (j*binwidth,i,j,-score,-1)  )
               j+=1
               while mask_iter_j<len(masked_segments) and masked_segments[mask_iter_j][1]<j*binwidth: mask_iter_j+=1

     print("#done making tile bumps.  len(tile_bumps)=",len(tile_bumps),strftime("%Y-%m-%d %H:%M:%S"),file=logfile,sep="\t")
     logfile.flush()
     if debug:
          for tile in tile_scores.keys():
               print("tile:",scaffold,tile[0],tile[1],tile[0]*binwidth,tile[1]*binwidth,tile_scores[tile],tile_counts[tile],model.tileScore(binwidth,tile,tile_counts[tile],tile_scores[tile],rangeCutoff=maxx,minSep=minx),sep="\t")
               for read in tile_reads[tile]:
                    print("tileread:",tile,read)

     tile_bumps.sort()
     print("#done sorting tile bumps",strftime("%Y-%m-%d %H:%M:%S"),file=logfile,sep="\t")
     logfile.flush()
     
     tripped=False
     row_sums={}
     col_sums={}
     row_counts={}
     col_counts={}
     break_buffer = []
     stretch_start=0
     minx=0
     state=0
     low_point=0
     ii=0
     last_trip=0
     while ii<len(tile_bumps):
          
          x,i,j,scoreD,dn = tile_bumps[ii]
          row_sums[i] = row_sums.get(i,0.0) + scoreD
          col_sums[j] = col_sums.get(j,0.0) + scoreD
          row_counts[i] = row_counts.get(i,0) + dn
          col_counts[j] = col_counts.get(j,0) + dn
          if dn==-1 and row_counts[i]==0: del row_sums[i] 
          if dn==-1 and col_counts[j]==0: del col_sums[j] 

          if ii%100000==0:
              print("#progress: ii= {} / {}".format(ii,len(tile_bumps)),x,i,j,scoreD,strftime("%Y-%m-%d %H:%M:%S"),file=logfile,sep="\t")
              logfile.flush()

          while ii<len(tile_bumps)-1 and tile_bumps[ii+1][0]==x:
               ii+=1
               x,i,j,scoreD,dn = tile_bumps[ii]
               row_sums[i] = row_sums.get(i,0.0) + scoreD
               col_sums[j] = col_sums.get(j,0.0) + scoreD
               row_counts[i] = row_counts.get(i,0) + dn
               col_counts[j] = col_counts.get(j,0) + dn
               if dn==-1 and row_counts[i]==0: del row_sums[i] 
               if dn==-1 and col_counts[j]==0: del col_sums[j] 

          total = sum(row_sums.values())

          row_vals = list(row_sums.values())
          col_vals = list(col_sums.values())

          row_vals.sort()
          col_vals.sort()

          without_best_row = sum(row_vals[:-nreject]) #total - max(row_sums.values())
          without_best_col = sum(col_vals[:-nreject]) #total - max(col_sums.values())
          trimmed_total = min(without_best_row,without_best_col,total)
          score=trimmed_total

          if (not clipped_hist==False):
              if min(x,slen-x)>scaffold_end_filter:
                  hist_bin=int(score/hist_bisize)*hist_bisize
                  clipped_hist[ hist_bin ] = clipped_hist.get(hist_bin,0)+1

          if score>t1:
               tripped=True
               last_trip=x
          if tripped and score<t2 and state==0:
               stretch_start=x
               state=1
               low_point =score
               minx = x
          if debug: print("trimmed_support",x,trimmed_total,total,without_best_row,without_best_col,tripped,maxtile*binwidth,scaffold,tripped,score<t2,state,low_point,minx,stretch_start,last_trip)
          if state==1 and score>t2:               
               if logfile: 
                    break_buffer.append( (scaffold,stretch_start,x,minx,low_point,slen)  )
#               print(scaffold,stretch_start,x,slen,low_point,"break")
               state=0
          if state==1:
               if score<low_point:
                    low_point=score
                    minx=x
          ii+=1
     print("#done building breaks buffer",strftime("%Y-%m-%d %H:%M:%S"),file=logfile,sep="\t")
     logfile.flush()

     if state==1:
         if logfile: 
             break_buffer.append( (scaffold,stretch_start,x,minx,low_point,slen)  )
         
     for scaffold,stretch_start,x,minx,low_point,slen in break_buffer:
          if (x < last_trip) or (minx < last_trip):
               min_fine_graph_segment_overlap_score=False
               for fg_scaffold,fg_stretch_start,fg_x,fg_minx,fg_low_point in fine_grain_segments:
                    if pairs_overlap((fg_stretch_start,fg_x),(stretch_start,x)):
                         if (not min_fine_graph_segment_overlap_score) or (min_fine_graph_segment_overlap_score > fg_low_point):
                              min_fine_graph_segment_overlap_score = fg_low_point
               logfile.write("{} {} {} {} {} {} {} clippedLLR\n".format(scaffold,stretch_start,x,minx,low_point,slen,min_fine_graph_segment_overlap_score))
     return fine_grain_support_curve
if __name__=="__main__":

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug',default=False,action="store_true",help="Turn on debugging ouput")
     parser.add_argument('-p','--progress',default=False,action="store_true",help="Print progress info")
     parser.add_argument('-l','--layout',default=False,help="File containing scaffolding layout.")
     parser.add_argument('-s','--segments',default=False,help="File containing scaffolding segments.")
#     parser.add_argument('-L','--length',default=False,help="File containing lenghts.")
     parser.add_argument('-m','--mask',default=False,help="File containing segments to mask.")
     parser.add_argument('-g','--plotfiles',default=False,help="Plot file name prefix.")

     parser.add_argument('-f','--logfile',default=False,help="Output file for storing score segments.")
     
     parser.add_argument('-q','--mapq',type=int,default=55)

     parser.add_argument('-t','--t1',type=float,default=50.0,help="Don't break in the trailing regions at the end of scaffolds where support never exceeds this number.")
     parser.add_argument('-T','--t2',type=float,default=25.0,help="Report all segments where support dips below this threshold, if they are not in the trailing ends.")

     parser.add_argument('-w','--binwidth',type=float,default=3000.0)
     parser.add_argument('-n','--nreject',type=int,default=2)
     parser.add_argument('-c','--my_chunk',type=int,default=1)
     parser.add_argument('-C','--nchunks',type=int,default=32)
     parser.add_argument('--minx',type=int,default=1000, help=" ")
     parser.add_argument('--maxx',type=int,default=200000,help=" ")
     parser.add_argument('-S','--scaffold',default=False)
     #     parser.add_argument('-b','--bamfiles',required=True)
     parser.add_argument('-b','--bamfile',action="append")
     parser.add_argument('-M','--model')

     args = parser.parse_args()
     if args.progress: print("#",args)

     fmodel=open( args.model )
     contents = fmodel.read()
     try:
          fit_params=eval(contents)
     except:
          print("couldn't set model parameters", args.model)
     fmodel.close
     ces.set_exp_insert_size_dist_fit_params(fit_params)

     if args.mask:
          mask_ranges = read_mask_ranges( open(args.mask) )
     else:
          mask_ranges={}
     
     #
     # Read in a hirise layout (from "^p:" lines)
     #
     mapper = MapperLite()
     mapper.load_layout(open(args.layout))

     if args.scaffold:
          my_scaffolds={args.scaffold:1}
     else:
          my_scaffolds={}
          scaffold_hashes={}
          for s in list(mapper.scaffolds.keys()):
               scaffold_hashes[s]=struct.unpack("<L", hashlib.md5(s.encode("utf-8")).digest()[:4])[0]%args.nchunks
               if scaffold_hashes[s]==args.my_chunk:
                    my_scaffolds[s]=1
                    if args.debug: print("my scaffold:",s)
          
     bamlist = [ pysam.Samfile(bamfile,"rb") for bamfile in args.bamfile ] 

     segment_logfile=False
     if args.logfile: segment_logfile = open(args.logfile,"wt")

     for sc in sorted(my_scaffolds.keys()):
          ncontigs=len(mapper.scaffold_contigs[sc]) 
          slen=mapper.scaffold_length[sc]
          if not ncontigs>1: continue
          print("sc:",sc,ncontigs,slen)
          #          print(sc,mapper.scaffold_contigs[sc])
          fp=False
          contigs = sorted(mapper.scaffold_contigs[sc],key=lambda x: mapper.contigx[x])
          ii=0
          gap_locations=[]
          for i in range(len(contigs)-1):
               con = contigs[i]
               x= max(mapper.contigx[con]+mapper.contig_strand[con]*mapper.contig_length[con],mapper.contigx[con])
               con2 = contigs[i+1]
               y= min(mapper.contigx[con2]+mapper.contig_strand[con2]*mapper.contig_length[con2],mapper.contigx[con2])
               gap_locations.append((int((x+y)/2),con,con2))

          if args.plotfiles: 
               fn="{}{}".format(args.plotfiles,sc)
               fp=open(fn,"wt")
               print("#plot \"{}\" i 2 u 1:2 w steps, \"\" i 1 u 1:($2/20) lt 3 pt 5 ps 0.7, \"\" i 0 u 1:(($2-5)*100) w steps lt 3, -500 lt 3".format(fn))
               ii=0
               fp.write("0\t0\n")
               for con in contigs:
                    ii+=1
                    ii=ii%2
                    if args.debug: print(con,mapper.contigx[con],mapper.contig_strand[con],mapper.contig_length[con],slen)
                    if mapper.contig_strand[con]==1:
                         fp.write("{}\t{}\n".format( mapper.contigx[con],2*ii-1 ))
                         fp.write("{}\t{}\n".format( mapper.contigx[con]+mapper.contig_length[con],0 ))
                    else:
                         fp.write("{}\t{}\n".format( mapper.contigx[con]-mapper.contig_length[con],2*ii-1 ))
                         fp.write("{}\t{}\n".format( mapper.contigx[con],0 ))
               fp.write("\n\n\n")

          masked_segments = make_scaffold_mask(sc,mapper,mask_ranges)
          pairs2support(sc,chicago_pairs(sc,mapper,bamlist,minq=args.mapq,mask=mask_ranges),ces.model,masked_segments=masked_segments,slen=mapper.scaffold_length[sc],buff=1,debug=args.debug,gpf=fp,joins=gap_locations,minx=args.minx,maxx=args.maxx,logfile=segment_logfile,t1=args.t1,t2=args.t2,binwidth=args.binwidth,nreject=args.nreject)
          if fp: fp.close()

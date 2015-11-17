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
import pysam
from bamtags import BamTags
import re

read_length=100

def mask_test(scaffold,x,mask_ranges):
     #print("#",scaffold,x,mask_ranges.get(scaffold))
     if scaffold in mask_ranges :
          for z,b in mask_ranges[ scaffold ]:
               if z<x and x + read_length < b:
                    #print("MASKED!")
                    return True
     return False

class SegMapper():
     def __init__(self,segment_file,segment_list=False):
          segment_info={}
          segments    ={}
          seg2len={}

          if segment_file:
               f=open(segment_file)
               while True:
                    l=f.readline()
                    if not l: break
                    dummy,scaffold,x,y,leng = l.strip().split()
                    x=int(x)
                    y=int(y)
                    leng=int(leng)
                    broken_name="{}_{}".format(scaffold,x)
                    segment_info[broken_name]=(scaffold,x,y)
                    seg2len[broken_name]=y-x
                    segments[scaffold] = segments.get(scaffold,[])+[(x,y,broken_name)]
               f.close()
          elif segment_list:
               for scaffold,x,y in segment_list:
                    broken_name="{}_{}".format(scaffold,x)
                    segment_info[broken_name]=(scaffold,x,y)
                    seg2len[broken_name]=y-x
                    segments[scaffold] = segments.get(scaffold,[])+[(x,y,broken_name)]                    
          else:
               raise Exception 

          for s in segments.keys():
               segments[s].sort()
          self.segment_info=segment_info
          self.segments=segments
          self.seg2len=seg2len

     def get_segment_info(self,seg):
#          if not seg in self.segment_info:
#               ((ocontig,base),) = re.findall("(.*)_(\d+)$",seg)
#               return(ocontig,base,1)
          scaffold,x,y=self.segment_info.get(seg,(seg,1,1))
          return (scaffold,x,y)

     def map_coord(self,c,yy):
          i=0
          if not c in self.segments:
               return (c+"_1",yy)
          while i<len(self.segments[c]) and yy>self.segments[c][i][1]: i+=1
          if i==len(self.segments[c]) and yy>self.segments[c][i-1][1]:
               print("couldn't find segment for",c,yy,self.segments[c])
               return(False,False)
#               raise Exception 
          z  = self.segments[c][i][0]
          c2 = self.segments[c][i][2]
          return(c2,1+yy-z)


def bam2chicagolinks(bamlist,my_contigs,smap,mask_ranges,mapq,minl,internal=False,tidlist=False):

#     print("bam2chicagolinks")
#     if not smap: smap=SegMapper("/dev/null")
#     print(my_contigs)
#     print(smap)
#     print(smap.segment_info)
     for seg in my_contigs.keys():
          if not smap:
               ti=bamlist[0].references.index(seg)
               l1 = bamlist[0].lengths[ti]
          links = {}
          tids = {}
          if smap:
               scaffold,x,y = smap.get_segment_info(seg)
               region = "{}:{}-{}".format(scaffold,x,y)
          else:
               region = seg
               x=-1
               y=-1
          print("#",seg,region)
          for bam in bamlist:
               #print(region)
               for a in bam.fetch(region=region): 
                    if a.mapq<mapq: continue
                    if a.is_duplicate: continue
                    if BamTags.mate_mapq(a)<mapq: continue
                    c1,c2,xx,yy = a.tid,a.rnext,a.pos,a.pnext
                    if xx<x: continue #print("#wtf?")
                    c1 = bam.getrname(c1)
                    c2 = bam.getrname(c2)

                    if mask_test(c1,xx,mask_ranges): continue
                    if mask_test(c2,yy,mask_ranges): continue
                    #print a

                    i=0
                    if smap:
                         c2,yyy = smap.map_coord(c2,yy)
                    else:
                         yyy=yy

                    if c2==seg and not internal : continue
                    if not smap and not c1==c2: continue
                    #print( seg,c2,1+xx-x,1+yy-z)
                    #if 1+xx-x<0: print ("#?",region,scaffold,1+xx-x,x,y,seg,c1,xx,x,yy,z)
                    links[c2] = links.get(c2,[]) + [( 1+xx-x,yyy  )]
                    tids[c2]  = tids.get(c2,[])  + [ a.query_name ]
          for c2 in links.keys():
               if False and not c2 in smap.seg2len:
                    print ("that's weird: {} not in table of segment lengths?".format(c2))
                    raise Exception 
               if smap and not c2 in smap.seg2len: continue
               if smap:
                    if (smap.seg2len.get(seg)>minl and smap.seg2len.get(c2)>minl): 
                         if tidlist:
                              yield( (seg,c2,smap.seg2len.get(seg),smap.seg2len.get(c2),len(links[c2]),links[c2], tids[c2]) )
                         else:
                              yield( (seg,c2,smap.seg2len.get(seg),smap.seg2len.get(c2),len(links[c2]),links[c2]) )
               else:
#                    print(seg,c2,l1,l1,len(links[c2]),links[c2])
                    if tidlist:
                         yield( (seg,c2,l1,l1,len(links[c2]),links[c2],tids[c2]) )
                    else:
                         yield( (seg,c2,l1,l1,len(links[c2]),links[c2]         ) )

def read_mask_ranges( f ):
     mask_ranges={}
     #f = open(args.mask)
     while True:
          l = f.readline()
          if not l: break
          if l[0]=="#": continue
          c=l.strip().split()
          mask_ranges[c[0]] = mask_ranges.get(c[0],[]) + [(int(c[1]),int(c[2]))] 
     f.close()
     for s in mask_ranges.keys():
          mask_ranges[s].sort()
     return mask_ranges

if __name__=="__main__":
     parser = argparse.ArgumentParser()

     #parser.add_argument('-i','--input')
     parser.add_argument('-d','--debug',default=False,action="store_true")
     parser.add_argument('-p','--progress',default=False,action="store_true")
     parser.add_argument('-L','--contiglist',default=False)
     parser.add_argument('-I','--internal',default=False,action="store_true")
     parser.add_argument('-S','--segments',default=False)
     parser.add_argument('-b','--bamfile',action="append")
     parser.add_argument('-R','--mask',default=False)
     parser.add_argument('-q','--mapq',default=10,type=int)
     parser.add_argument('-m','--minl',default=1000,type=int)


     args = parser.parse_args()
     if args.progress: print("#",args)

     if args.segments:
          smap=SegMapper(args.segments)
     else:
          smap=False

     if args.mask:
          mask_ranges = read_mask_ranges( open(args.mask) )
     else:
          mask_ranges={}

     my_contigs={}
     if args.contiglist:
          f=open(args.contiglist)
          while True:
               l=f.readline()
               if not l: break
               contig = l.strip().split()[0]
               my_contigs[contig]=1
          f.close()

     bamlist = [ pysam.Samfile(bamfile,"rb") for bamfile in args.bamfile ] 
     for links_tuple in bam2chicagolinks(bamlist,my_contigs,smap,mask_ranges,args.mapq,args.minl,args.internal):
          print("\t".join(map(str,links_tuple)))


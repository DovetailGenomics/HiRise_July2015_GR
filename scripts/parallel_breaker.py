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
from hirise_assembly import HiriseAssembly
import multiprocessing
from multiprocessing import Pool, Queue, Process, JoinableQueue
import chicago_edge_scores as ces
from chicago_support_bootstrap import pairs2support
import gc
import os
from time import gmtime, strftime,sleep
#import psutil

def pair2bin2D(pair,binwidth):
     a,b=pair
     return( int(a/binwidth),int(b/binwidth) )

#def pairs2support(scaffold,pairs,model,slen=0,masked_segments=[],mapper=None,buff=500,debug=False,t1=20.0,t2=5.0,gpf=False,minx=1000,maxx=1e7,joins=[],logfile=False,binwidth=1000,nreject=2):

def handle_scaffold(q,histogram_queue,myid,outfile,model,hra,t1=5,t2=20,debug=False,apply_mask=False,maxx=200000,minx=500,binwidth=1500,nreject=2):
     ns=[0,0,0]
     buff=0
     while True:
          raw_support_histogram={}
          clipped_support_histogram={}
          scaffold,slen,pairs = q.get()
#          print("#")
          outfile.write("#{}\t{}\t{}\n".format(scaffold,slen,len(pairs)))
          mask=[]
          if apply_mask:
               mask = hra.make_scaffold_mask(scaffold)
          pairs2support(scaffold,pairs,model,slen,logfile=outfile,masked_segments=mask,t1=t1,t2=t2,binwidth=binwidth,maxx=maxx,minx=minx,buff=0,nreject=nreject,debug=debug,raw_hist=raw_support_histogram, clipped_hist=clipped_support_histogram )
          outfile.flush()
          q.task_done()
          outfile.write("#{} done\n".format(scaffold))
          outfile.flush()
          histogram_queue.put((raw_support_histogram,clipped_support_histogram))

def merge_histograms(q,final_h):
     raw_support_histogram={}
     clipped_support_histogram={}
     print("# merge histogram worker started")
     while True:
          (raw_support_histogram_delta,clipped_support_histogram_delta) = q.get()
          if raw_support_histogram_delta==-999:
               print("#got histogram merge finalize signal")
               sys.stdout.flush()
               final_h.put((raw_support_histogram,clipped_support_histogram))
               final_h.close()
               q.task_done()
               break
          else:
               #print("#merge:",len(raw_support_histogram_delta.keys()))
               sys.stdout.flush()
               for b,c in raw_support_histogram_delta.items():
                    raw_support_histogram[b] = raw_support_histogram.get(b,0)+c
               for b,c in clipped_support_histogram_delta.items():
                    clipped_support_histogram[b] = clipped_support_histogram.get(b,0)+c
               q.task_done()

def read_pairs_from_bam(inq,i,bamfile,hra,mapq,scaffold_slice=False):
     #print("x",i,bamfile)
     ns=[0,0,0]
     buff=0

     pair_buffer={}

     scaffold_contig_counts={}
     contig_scaffold_counts={}

     scaffolds = list(hra.scaffold_lengths.keys())
     scaffolds.sort(key=lambda x: (hra.scaffold_lengths[x],x),reverse=True)
     my_scaffolds=False
     if not scaffold_slice==False:
          my_scaffolds=[]
          slicei,slicen = scaffold_slice.split(",")
          slicei,slicen = int(slicei),int(slicen)
          for i,scaffold in enumerate(scaffolds):
               if i%slicen == slicei:
                    my_scaffolds.append(scaffold)

     for c,s in hra.contig_scaffold.items():
          contig,base = c
          if (not my_scaffolds==False) and (not s in my_scaffolds): continue 
          if not s in scaffold_contig_counts: 
               scaffold_contig_counts[s]={}
          if not contig in contig_scaffold_counts:
               contig_scaffold_counts[contig]={}

          scaffold_contig_counts[s][contig] = scaffold_contig_counts[s].get(contig,0)+1
          contig_scaffold_counts[contig][s] = contig_scaffold_counts[contig].get(s,0)+1

          #print(s,contig,scaffold_contig_counts[s][contig])

#     def process_scaffold(scaffold):
#          print("process:",scaffold,len(pair_buffer.get(scaffold,[])),pair_buffer.get(scaffold,[]),sep="\t")
#          print("process:",scaffold,len(pair_buffer.get(scaffold,[])),sep="\t")
#          q.put((scaffold,hra.scaffold_lengths[scaffold],pair_buffer.get(scaffold,[])))
#          print(multiprocessing.active_children())
#          if scaffold in pair_buffer: del pair_buffer[scaffold]

#     gc.disable()
     for row in hra.chicago_pairs(scaffold_contig_counts=scaffold_contig_counts,contig_scaffold_counts=contig_scaffold_counts,callback=False,bamfile=bamfile,mapq=mapq,my_scaffolds=my_scaffolds):
          if len(row)==7:
               s1,s2,a,b,c1,c2,tid=row
          else:
               #print("done with scaffold:",i,row)
               done_scaffold=row[0]
               inq.put( (i,done_scaffold,pair_buffer.pop(done_scaffold,[])  ) )
               continue
          if not s1==s2: continue
          
#          print(s1,s2,x,y,scaffold_contig_counts.get(s1),sep="\t")
          if not s1 in pair_buffer: 
               pair_buffer[s1]=[]
#          pair_buffer[s1].append((s1,a  ,b  ,c1,c2,0,  0,     0      , s1,     0,   0,  0,   s2,     0,   0, 0  ,tid) )
          pair_buffer[s1].append((a,b))

     remaining_scaffolds = list(pair_buffer.keys())
     for done_scaffold in remaining_scaffolds:
          print("trailing scaffold:",i,done_scaffold,len(pair_buffer[done_scaffold]))
          inq.put( (i,done_scaffold,pair_buffer.pop(done_scaffold,[])  ) )

     inq.put( (i,"DONE",[] ) )
     inq.close()
     print("#reader done",i,bamfile)
     return(0)

if __name__=="__main__":

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug',default=False,action="store_true",help="Turn on debugging ouput")
     parser.add_argument('-p','--progress',default=False,action="store_true",help="Print progress info")
     parser.add_argument('-M','--nomask',default=False,action="store_true",help="Don't apply masking correction.")
     parser.add_argument('--top',default=False,action="store_true",help="Report cpu loads")

     parser.add_argument('-t','--t1',type=float,default=5.0,help="Don't break in the trailing regions at the end of scaffolds where support never exceeds this number.")
     parser.add_argument('-T','--t2',type=float,default=20.0,help="Report all segments where support dips below this threshold, if they are not in the trailing ends.")
     parser.add_argument('-q','--mapq',type=int,default=55)

     parser.add_argument('-S','--slice',default=False)

     parser.add_argument('-i','--infile',help="File containing .hra formatted hirise assembly.")
     parser.add_argument('-o','--outfile',help="Output file name.")
     parser.add_argument('-H','--histogram_file',help="Save score histograms to FILE")
     parser.add_argument('-j','--nthreads',help="Number of threads.",type=int,default=16)
     

     args = parser.parse_args()

     if args.top:
          import psutil

     hra = HiriseAssembly()
     hra.load_assembly(args.infile)

     hra.merge_masked_regions(debug=args.debug)

     print(len(hra.layout_lines))
     if len(hra.layout_lines)==0: 
          print("#make trivial layout")
          hra.make_trivial_layout_lines(debug=args.debug)

     ces.set_exp_insert_size_dist_fit_params(hra.model_params)
     model=ces.model

     nbams = len(hra.bams)
     inq = JoinableQueue(maxsize=0)
     readers=[]
     for i in range(nbams):
          reader = Process(target=read_pairs_from_bam, args=(inq,i,hra.bams[i],hra,args.mapq,args.slice),daemon=False)
          reader.start()
          readers.append(reader)

     q = JoinableQueue(maxsize=0)
     histogram_queue = JoinableQueue(maxsize=0)
     final_histogram_queue = JoinableQueue(maxsize=0)
     outfiles = [ open("{}.part.{}".format(args.outfile,i),"wt") for i in range(args.nthreads) ]

     workers=[]
     for i in range(args.nthreads):
          worker = Process(target=handle_scaffold, args=(q,histogram_queue,i,outfiles[i],model,hra,args.t1,args.t2,args.debug,not args.nomask),daemon=True)
          worker.start()
          workers.append(worker)

     histogram_merge_worker = Process(target=merge_histograms, args=(histogram_queue,final_histogram_queue),daemon=True)
     histogram_merge_worker.start()

     if args.top:
          reader_procs = [ psutil.Process(reader.pid) for reader in readers ]
          worker_procs = [ psutil.Process(worker.pid) for worker in workers ]

     pair_buffer={}
     scaffold_count={}
#     while (not inq.empty()) or sum( [reader.is_alive() for reader in readers] )>0:
     while True:
          if args.debug: print("get")
          try:
               procid,scaffold,pairs = inq.get()
#               procid,scaffold,pairs = inq.get(True,10)
               #print("#got data:",procid,scaffold,len(pairs))
               print("#got data from inq:",procid,scaffold,len(pairs),inq.empty(),inq.qsize(),inq.full(),strftime("%Y-%m-%d %H:%M:%S"),sum( [reader.is_alive() for reader in readers] ),"q.size():",q.qsize(),file=sys.stderr,sep="\t")
               sys.stderr.flush()
               sys.stdout.flush()
          except Exception as e:
               print(e,file=sys.stderr)
               if args.top:
                    print("queue get timed out",[reader.cpu_percent() for reader in reader_procs],[worker.cpu_percent() for worker in worker_procs])
               #print("#timed out",inq.empty())
               print("#read from queue timed out:",inq.empty(),inq.qsize(),inq.full(),strftime("%Y-%m-%d %H:%M:%S"),sum( [reader.is_alive() for reader in readers] ),file=sys.stderr,sep="\t")
               sys.stderr.flush()
               continue
          if args.debug: print("got")
          if not scaffold in pair_buffer:
               pair_buffer[scaffold]=[]
          pair_buffer[scaffold] += pairs
          scaffold_count[scaffold] = scaffold_count.get(scaffold,0)+1
#          print("#:",scaffold,scaffold_count[scaffold])
          if scaffold=="DONE": 
               print("#seen DONE {} times.".format(scaffold_count.get(scaffold,0)),strftime("%Y-%m-%d %H:%M:%S"),file=sys.stderr,sep="\t")
               if scaffold_count[scaffold] == nbams:
                    print("#done. calling task_done:",strftime("%Y-%m-%d %H:%M:%S"),file=sys.stderr,sep="\t")
                    sys.stderr.flush()
                    inq.task_done()
                    print("#done. called task_done(); break.",strftime("%Y-%m-%d %H:%M:%S"),file=sys.stderr,sep="\t")
                    sys.stderr.flush()
                    break
               else:
                    print("#x: not enought {}".format(scaffold_count[scaffold]),file=sys.stderr,sep="\t")
                    inq.task_done()
          else:
               if args.debug: print("##############3",procid,scaffold,scaffold_count.get(scaffold),len(pairs),len(pair_buffer[scaffold]))

               if scaffold_count[scaffold] == nbams:
     #               print("#q.put ... ",scaffold,len(pair_buffer.get(scaffold,[])))
     #               if args.debug:
     #                    for pp in pair_buffer[scaffold]:
     #                         print("xy:",scaffold,min(pp[1],pp[2]),max(pp[1],pp[2]))
                    q.put((scaffold,hra.scaffold_lengths[scaffold],pair_buffer.pop(scaffold,[])))
     #               print("#.")

     #          pairs2support(scaffold,pairs,model,slen,logfile=outfile)
     #          outfile.flush()
               if args.debug: print([reader.is_alive() for reader in readers])
     #          if sum( [reader.is_alive() for reader in readers ] )==0: break
               inq.task_done()

     print("#Outside main loop.",strftime("%Y-%m-%d %H:%M:%S"),file=sys.stderr,sep="\t")
     print("#inq empty?",inq.empty())
     print("#waiting for worker threads to complete. ({})".format( sum([ worker.is_alive() for worker in workers ]) ) )
     sys.stdout.flush()
     sys.stderr.flush()
     while True:
          print("#Scaffold in the queue =",q.qsize(), "workers running:", sum([ worker.is_alive() for worker in workers ]) ,strftime("%Y-%m-%d %H:%M:%S"),histogram_queue.qsize(),file=sys.stderr,sep="\t")
          sys.stderr.flush()
          sleep(10)
          if q.qsize()==0: break

     q.join()
     print("#workers done.")
     sys.stdout.flush()
     print("#All worker threads joined.",strftime("%Y-%m-%d %H:%M:%S"),file=sys.stderr,sep="\t")
     sys.stderr.flush()
     
#     inq.join()
     for fh in outfiles:
#          print(fh)
          fh.flush()
          fh.close()
     
     of=open(args.outfile,"wt")
     for outfile in [ "{}.part.{}".format(args.outfile,i) for i in range(args.nthreads) ]:
          for l in open(outfile,"rt"):
               of.write(l)
          os.remove(outfile)
     of.close()
     
     print("#Finalize histogram merging.",strftime("%Y-%m-%d %H:%M:%S"),file=sys.stderr,sep="\t")
     sys.stderr.flush()
     histogram_queue.put((-999,0))
     histogram_queue.close()

     print("#Get merged histograms.",strftime("%Y-%m-%d %H:%M:%S"),file=sys.stderr,sep="\t")
     sys.stderr.flush()
     (raw_support_histogram,clipped_support_histogram) = final_histogram_queue.get()

     if args.histogram_file:
          hf=open(args.histogram_file,"wt")
          for b,c in sorted(raw_support_histogram.items()):
               print("raw:",b,c,sep="\t",file=hf)
          for b,c in sorted(clipped_support_histogram.items()):
               print("clipped:",b,c,sep="\t",file=hf)
          hf.close()

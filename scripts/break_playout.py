#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
from hirise_assembly import HiriseAssembly
import shlex
import shutil
import subprocess
import tempfile

def make_break_fig(x):
     scaffold,a,b,c,score,td,mapq,maskfile = x
     if not maskfile: maskfile="/dev/null"
     ticks = tempfile.NamedTemporaryFile(dir=td,suffix=".txt",delete=False)
     tf    = tempfile.NamedTemporaryFile(dir=td,suffix=".pdf",delete=False)
     ticks.write("{scaffold}\t{a}\n".format(scaffold=scaffold,a=c).encode())
     ticks.flush()
     cmd="oovis.py -T {tickfile} -i {layout} -R {scaffold}:{a}-{b} -p 300000 --noreference -o {tf} -f 1000 -d --rule 50000 -D -I -b {bams} -q {mapq} --mask {maskfile}".format(mapq=int(mapq), bams=" -b ".join(args.bams), layout=args.infile, scaffold=scaffold, a=a, b=b, tf=tf.name , tickfile=ticks.name, maskfile=maskfile)
     print(cmd)
     #if args.debug: print(cmd)
     output,error = subprocess.Popen(shlex.split(cmd),stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate()
     if True or args.debug:
          print(output)
          print(error)

     tf.flush()
     tf.close() 
     header="{scaffold}:{a} {score}".format(scaffold=scaffold,a=c,score=score)
     return (tf.name,header)

def testf(x):
     a,b,c=x
     print(a,b,c)

if __name__=="__main__":
     import sys
     import argparse

     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     parser.add_argument('-d','--debug',default=False  ,action="store_true",help="Turn on debugging ouput")
     parser.add_argument('-b','--breaks',default=False ,help="File containing breaks")
     parser.add_argument('-B','--bams',default=[],action="append" ,help="File containing breaks")
     parser.add_argument('-i','--infile',default=False ,help="Input layout in p: lines format.")
     parser.add_argument('-j','--nprocesses',default=16,type=int ,help="Number of subprocesses to use when making PDF reports.")
     parser.add_argument('-H','--head',default=False,type=int ,help="Maximum number of breaks to make visualizations for.")
     parser.add_argument('-o','--outfile',default=False,help="Filename for output.")
     parser.add_argument('-R','--report',default=False ,metavar="FILE", help="Generate a PDF report with visualizations of the breaks made.")
     parser.add_argument('-t','--threshold',default=5.0,type=float,help="Minimum acceptable LLR support.")
     parser.add_argument('-T','--contig_break_threshold',type=float,help="Minimum acceptable LLR support within contig sequence.")
     parser.add_argument('-q','--mapq',default=10.0,type=float,help="Mapping quality score cutoff.")
     parser.add_argument('--mask',default=False,help="File with segments to mask for dot-plot visualization.")

     args = parser.parse_args()
     if args.contig_break_threshold==None: args.contig_break_threshold = args.threshold

     asf = HiriseAssembly()
     asf.load_playout(args.infile)

     breaks=[]
     scores={}
     stype={}
     for l in open(args.breaks):
          if l[0]=="#": continue
          fields = l.strip().split()
          scaffold,a,b,c,score,slen = fields[:6]

          lowest_raw_score=args.threshold+5
          clipped=False
          if "clippedLLR" in fields:
               clipped=True
               if not fields[6]=="False":
                    lowest_raw_score = float(fields[6])
#10042 241876 242439 241876 13.357019490186218 341715 rawLLR
#10042 243000.0 246000.0 243000.0 11.97365670913224 341715 False clippedLLR

          if clipped and lowest_raw_score < args.threshold: continue # this is redundant with a fine-grained break we've alreday got.

          a=int(float(a)) #start
          b=int(float(b)) #end
          c=int(float(c)) #lowpoint
          score=float(score)
          #if b-a > 1000:
          #     breaks.append((scaffold,a))
          #     breaks.append((scaffold,b))
          #else:
          if score<args.threshold:
               breaks.append((scaffold,a,b,c))
               scores[scaffold,a,b,c]=score
               if clipped: 
                    stype[scaffold,a,b,c]="clipped"

     if args.report:
          import re
          from multiprocessing import Pool
          from functools import partial
          import random

          td=tempfile.mkdtemp()
          tfs=[]
          headers=[]

          p=Pool(args.nprocesses)
          if args.head:
               random.shuffle(breaks)
               tfs_headers=p.map( make_break_fig, [ (scaffold,a,b,c,scores[scaffold,a,b,c],td,args.mapq,args.mask ) for scaffold,a,b,c, in breaks[:args.head]])
          else:
               tfs_headers=p.map( make_break_fig, [ (scaffold,a,b,c,scores[scaffold,a,b,c],td,args.mapq,args.mask ) for scaffold,a,b,c, in breaks ] )


          print(tfs_headers)
          orgfile=tempfile.NamedTemporaryFile(dir=td,suffix=".org")
          orgfile.write("""#+OPTIONS:   H:0 num:t toc:nil \\n:nil @:t ::t |:t ^:nil -:t f:t *:t <:t
#+LaTeX_CLASS: article
#+LaTeX_CLASS_OPTIONS: [legalpaper,landscape,utopia,8pt,listings-sv,microtype,paralist,DIV=15,BCOR=15mm]
#+LaTeX_HEADER: \\usepackage[legalpaper, total={{12in, 7in}}]{{geometry}}

#+TITLE: Breaks on {layout} for {title}

""".format(layout=args.infile,title=args.breaks).encode())

          #print(tfs)

          for tf,header in tfs_headers:
               orgfile.write("""* Break {label}

\\includegraphics[width=70em]{{{pdffile}}}


""".format(pdffile=tf,label=header).encode())
          orgfile.flush()
          
          cmd = "emacs {} --batch -f org-export-as-pdf  --kill ".format(orgfile.name)
          print(cmd)
          output,error = subprocess.Popen(shlex.split(cmd),stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate()
          if args.debug:
               print(output)
               print(error)

          of=re.sub(".org",".pdf",orgfile.name)
          cmd = "mv {} {}".format(of,args.report)
          output,error = subprocess.Popen(shlex.split(cmd),stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate()
          if args.debug:
               print(output)
               print(error)
          

     if args.outfile:
          asf.add_breakpoints_ranges(breaks,debug=args.debug,scores=scores,unlink_threshold=args.threshold, contig_break_threshold=args.contig_break_threshold,score_types=stype)
          asf.validate()
          asf.dump_layout( open(args.outfile,"wt") )


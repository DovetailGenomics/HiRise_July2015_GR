#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from past.utils import old_div
from builtins import object
import sys
import argparse
import networkx as nx
import sys

parser = argparse.ArgumentParser()

#parser.add_argument('-i','--input')
parser.add_argument('-d','--debug',default=False,action="store_true")
parser.add_argument('-B','--bridges',default=False,action="store_true")
parser.add_argument('-F','--fakeedges',default=False,action="store_true")
parser.add_argument('-b','--blacklist')
parser.add_argument('-L','--links',default=[],action="append")
parser.add_argument('-t','--t1',default=False,type=float)
parser.add_argument('-D','--clusterMap',default=False)
parser.add_argument('-v','--t2',default=25.0,type=float,help="threshold used in defining and removing bridges")
parser.add_argument('-p','--progress',default=False,action="store_true")

args = parser.parse_args()
if args.progress: print("#",args)

g=nx.Graph()

rl = sys.setrecursionlimit(50000)

bl={}
if args.blacklist:
  f=open(args.blacklist)
  while True:
    l=f.readline()
    if not l: break
    c=l.strip().split()
    bl[c[0]]=1

def bridges(G):
  """
  Bridge detection algorithm based on WGB-DFS.

  """

  if G.is_directed():
    raise nx.NetworkXError('This function is for undirected graphs.\n'
                           'Use directed_wgb_dfs() for directed graphs.')

  class WhiteGrayBlackDFS(object):
    def __init__(self, G):
      # white: empty
      # gray: 1
      # black: 2

      self.visited = set()
      self.dfs_num = {}
      self.num = 0
      self.G = G
      self.back_edges = {} #defaultdict(set)

    def bridges(self, parent, current):
      #print '{'
      #print 'parent, current:', parent, current
      #print 'dfs_num:', self.dfs_num
      self.visited.add(current)
      current_lowpoint = self.dfs_num[current] = self.num

      self.num += 1
      #print 'dfs_num:', self.dfs_num

      for child in G.neighbors(current):
        if child != parent:
          #print 'current, child:', current, child
          if not current in self.back_edges or (current in self.back_edges and not child in self.back_edges[current]):
            if child in self.visited:
              current_lowpoint = min(current_lowpoint, self.dfs_num[child])
            else:
              for x in self.bridges(current, child):
                yield x
              if self.child_lowpoint > self.dfs_num[current]:
                #print '>>> bridge:', current, child
                yield (current, child)
              current_lowpoint = min(current_lowpoint, self.child_lowpoint)

      #print 'parent, current, current_lowpoint:', parent, current, current_lowpoint
      #print 'dfs_num:', self.dfs_num
      #print '}'
      self.child_lowpoint = current_lowpoint


  dfs = WhiteGrayBlackDFS(G)

  for x in G:
    if not x in dfs.visited: 
      #print x
      for e in dfs.bridges(x, x):
        yield e

weight={}                                   
while True:
     l=sys.stdin.readline()
     if not l: break
     if l[0]=="#": continue
     c=l.strip().split()
     if args.t1 and float(c[2])<args.t1: continue
     if (c[0] in bl) or (c[1] in bl): continue
     g.add_edge(c[0],c[1])
     weight[c[0],c[1]]=float(c[2])
     weight[c[1],c[0]]=float(c[2])
nb=0
non_cherry_bridges=0

if args.bridges:
  to_remove=[]
  for e in bridges(g):
       #print e
       if (not g.degree(e[0])==1) and (not g.degree(e[1])==1): 
            non_cherry_bridges+=1
            if weight[e]<args.t2: 
                 to_remove.append(e)
       nb+=1

  ne=0
  for e in g.edges():
       ne+=1

  print("bridges:",non_cherry_bridges,nb,ne)


  g.remove_edges_from(to_remove)

fd=False
if args.clusterMap:
  fd=open(args.clusterMap,"w")

i=0
lll=[]
clustermap={}
clustersize={}

for c in nx.connected_components(g):
  i+=1
  clustersize[i]=len(c)
#  lll.append(len(c))
  if fd:
    for cc in c:
      fd.write("{}\t{}\n".format(cc,i))
  for cc in c: clustermap[cc]=i
if fd: fd.close()

import glob
node2clustercounts={}
for linkfile_glob in args.links:
  for f in glob.glob(linkfile_glob):
    print("#",f)
    fh=open(f,"rt")
    while True:
      l=fh.readline()
      if not l: break
      if l[0]=="#": continue
      c=l.strip().split()
      if not( int(c[2])>1000 and int(c[3])>1000): continue

      if clustersize.get(clustermap.get(c[0]),1)==1:
        if not c[0] in node2clustercounts: node2clustercounts[c[0]]={}
        if c[1] in clustermap:
          node2clustercounts[c[0]][clustermap[c[1]]] = node2clustercounts[c[0]].get(clustermap[c[1]],0)+int(c[4])

for singleton in node2clustercounts.keys():
    others = list(node2clustercounts[singleton].keys())
    others.sort(key=lambda x:node2clustercounts[singleton][x],reverse=True)
    candidates =  [(x,node2clustercounts[singleton][x]) for x in others if node2clustercounts[singleton][x]>1]

    if len(candidates)==0: continue
    if len(candidates)==1 or candidates[0][1]>2*candidates[1][1]:
      clustermap[singleton]=candidates[0][0]
      clustersize[candidates[0][0]]+=1

lll=list(clustersize.values())
lll.sort(reverse=True)

n=sum(lll)
if len(lll)>1:
  print("#",n,len(lll),old_div(float(lll[0]),n),old_div(float(lll[0]),lll[1]),lll[:10])
else:
  print("#",n,len(lll),old_div(float(lll[0]),n),"nan",lll)
  
if args.fakeedges:
  exemplar={}
  for n,c in clustermap.items():
    if not c in exemplar:
      exemplar[c]=n
    else:
      print("{}\t{}\t{}".format(n,exemplar[c],1))

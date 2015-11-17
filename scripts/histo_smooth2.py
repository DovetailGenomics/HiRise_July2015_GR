#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import map
from past.utils import old_div
import sys

sc = 1000.0
if len(sys.argv)>1:
    sc=float(sys.argv[1])

queue=[]

bin_start = False
bin_fill = 0.0
ys=0.0
xs=0.0
nb=0.0
lastx=False
while True:
    l = sys.stdin.readline()
    if not l: break
    c = list(map(float,l.strip().split()))
    if not lastx:
        lastx=c[0]-1
    ys += c[1]
    xs += c[0]*c[1]
    nb += c[0]-lastx
    if ys > sc:
        print(old_div(xs,ys), old_div(ys,nb), nb)
        ys = 0.0
        xs = 0.0
        nb = 0.0
    lastx=c[0]

print(old_div(xs,ys), old_div(ys,nb), nb)

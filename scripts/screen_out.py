#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
import sys
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('rejectlist',metavar='REJECT',help="File with a list of things to filter out.")
parser.add_argument('column',type=int,metavar='COLUMN',help="Column of stdin to filter.")
parser.add_argument('-k','--keep',action="store_true",default=False,help="Keep them instead of filtering them out.")

args = parser.parse_args()

#print(args)

f=open(args.rejectlist)

keep=args.keep

reject_list={}
while True:
    l=f.readline()
    if not l:
        break
    reject_list[l.strip()]=1
f.close()

ci = args.column-1
#print "#",ci
while True:
    l=sys.stdin.readline()

    if not l:
        break

    if l[0]=="#": continue
    c=l.strip().split()
    if ((not keep) and (c[ci] not in reject_list)) or (keep and c[ci] in reject_list):
        print(l.strip())

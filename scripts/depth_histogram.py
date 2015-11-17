#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#


## Double threshold

import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-o','--outhist')
parser.add_argument("-t", "--threshold")
args = parser.parse_args()

h={}

while True:
    l= sys.stdin.readline()
    if not l: break
    c = l.strip().split()
    scaffold,x,y = c[0],int(c[1]),int(c[2])
    h[y]=h.get(y,0)+1

depths=list(h.keys())
depths.sort()
total_sum = 0;
for d in depths:    
    total_sum = total_sum + h[d];

if args.outhist:
    fh=open(args.outhist,"wt")
    for d in depths:
        fh.write("{}\t{}\n".format(d,h[d]))
    fh.close()

if args.threshold:
    fh = open(args.threshold, "wt")
    err = open("depth.err", "wt")
    depths=list(h.keys())
    depths.sort()
    cum_sum = [];
    last = 0; 
    q25 = -1;
    q75 = -1;
    # Loop through and calculate the cumulative sum of the data
    # recording the 25th and 75th quantiles as we find them.
    for d in depths :
        current = last + h[d];
        quantile = (1.0 * current)/total_sum;
        cum_sum.append(quantile);
        if q25 < 0 and quantile >= .25 :
            q25 = d;
        if q75 < 0 and quantile >= .75 :
            q75 = d;
        last = current;
    
    # Set threshold based on the inter quantile range (IQR) and
    # 75th threshold. So our outlier threshold is based up on the spread
    # and uppwer quantile of the data.
    # John Tukey, Exploratory Data Analysis, Addison-Wesley, 1977, pp. 43-44. 
    iqr = q75 - q25;
    print("Q25 is %d and Q75 is %d" % (q25,q75));
    fh.write("{}\n".format(q75 + 3 * iqr));
    fh.write("{}\n".format(q75 + 3.5 * iqr));
    fh.close()
    err.close()

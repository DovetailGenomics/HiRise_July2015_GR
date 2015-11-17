#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#


from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import sys
import string
import random
tt = str.maketrans("ACTGactg","TGACtgac")

class fastaWriter(object):
    def __init__(self,filename,linelen=60):
#        print(filename, file=sys.stderr)
        self.f=open(filename,"w")
#        print(self.f, file=sys.stderr)
        self.linelen=linelen
        self.buff=""
        self.x=0
        self.name=""

    def write(self,s):
        self.buff+=s
        self.x+=len(s)
#        print(len(self.buff),self.linelen, file=sys.stderr)
        #        while len(self.buff)>self.linelen:
        wrtn=0
        for i in range( 0, self.linelen*(int(old_div(len(self.buff),self.linelen))) , self.linelen ):
            wrtn+=self.linelen
            self.f.write( self.buff[i:i+self.linelen]+"\n" )
#                print "#",self.buff[i:i+self.linelen]
#                sys.stdout.flush()
        if (len(self.buff)%self.linelen)==0:
            self.buff=""
        else:
            self.buff= self.buff[ -(len(self.buff)%self.linelen) :]
#        print(len(self.buff),wrtn, file=sys.stderr)
        self.f.flush()
        
    def flush(self):
#        while len(self.buff)>self.linelen:
#        print(len(self.buff),self.linelen,"flush", file=sys.stderr)
        wrtn=0
        for i in range( 0, self.linelen*(int(old_div(len(self.buff),self.linelen))) , self.linelen ):
            wrtn+=self.linelen
            self.f.write( self.buff[i:i+self.linelen]+"\n" )
        if len(self.buff)%self.linelen>0:
            wrtn+=(len(self.buff)%self.linelen)
#            print("flush",self.buff[ -(len(self.buff)%self.linelen) :], file=sys.stderr)
            self.f.write( self.buff[ -(len(self.buff)%self.linelen) :] +"\n")
        self.buff=""
#        print(len(self.buff),wrtn, file=sys.stderr)
        self.f.flush()
        
    def next(self,name):
#        if self.x>0:
            
            #            print("#.",self.name,self.x, file=sys.stderr)
        self.x=0
        self.flush()
#        print(name, file=sys.stderr)
        sys.stdout.flush()
        self.f.write(">{}\n".format(name))
        self.f.flush()
        self.name=name

    def close(self):
        self.f.close()

def rc(s):
    s = s.translate(tt)
    return s[::-1]

def slurp_fasta(x):
    seqs={}
    f=open(x)
    b=""
    name=""
    while True:
        l=f.readline()
        if not l: break
        if l[0]==">":
            if len(b)>0:
                seqs[name]=b
                b=""
            c=l[1:].strip().split()
            name=c[0]
        else:
            b=b+l.strip()
    seqs[name]=b
    f.close()
    return seqs

if __name__=="__main__":
    import sys

#    print " ".join(sys.argv)

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-f','--fasta',required=True)
    parser.add_argument('-o','--outfile',required=True)
    parser.add_argument('-p','--prefix',required=False)
    parser.add_argument('-g','--gapclose',required=False)
    parser.add_argument('-d','--debug',default=False,action="store_true")
    parser.add_argument('--seed',required=False,type=int,default=1, help="Seed for random number generation, use -1 for no seed")
    args = parser.parse_args()

    if args.seed != -1 :
        random.seed(args.seed)

    fastafile=args.fasta
    seqs=slurp_fasta(fastafile)
#    print("slurped fastafile", file=sys.stderr)

    outfasta=fastaWriter(args.outfile)
#    print("opened outfile", file=sys.stderr)

    name_prefix=""
    if args.prefix:
        name_prefix=args.prefix
    else:
        import idGen as idGen
        name_prefix="Sc" + idGen.id()
        print("name_prefix",name_prefix,file=sys.stderr)
    #    for i in range(5):
    #        name_prefix += random.choice( "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz" )

#31      Scaffold142609_1.5      GGGGTACCAGGTGTCTGGGTGCTGGCAGAGCCGGCACTGATGTTTTCTGGG     Scaffold41854_1.5       GGTACCAGGTGTCTGGGTGCTGGCAGAGCCGGCACTGATGTTTTCTGGGGG     GGGGTACCAGGTGTCTGGGTGCTGGCAGAGCCGGCACTGATGTTTTCTGGGGG

    gapclose_data={}
    if args.gapclose:
        fh=open(args.gapclose)
        while True:
            l=fh.readline()
            if not l: break
            if l[0]=="#": continue
            scaffold,end1,kmer1,end2,kmer2,closure=l.strip().split()
            gapclose_data[end1,end2]=(kmer1,kmer2,closure)
        fh.close()

    last_patch=""
    kmer=51
    trailing_kemr=""
    nibbled_seq=""
    last_c=-1
    last_x=0
    last_s=0
    cl=0
    while True:
        l=sys.stdin.readline()
        if not l: break
        if not l[:2]=="p:": continue
#        print(l.strip(), file=sys.stderr)
        c=l.strip().split()
        #p: 1 Scaffold63081_1.5 25308872.0 8 31235809 37446086.0 ['79538.0', '8', '+', '31235809', '31278771', '1', '42972', '42972']
        scaffold,end,x = int(c[1]),c[2],float(c[3])
        if scaffold>last_s:
            outfasta.next(name_prefix+"_"+str(scaffold))
            last_x=0
        contig,p=end[:-2],end[-1:]
#        print(scaffold,end,x,contig,p, file=sys.stderr)
        if not contig == last_c:

            nibble_off=0
            nibbled_seq = ""
            if last_x>0:
                if (last_end,end) in gapclose_data:
                    kmer1,kmer2,closure= gapclose_data[last_end,end]
                    if not kmer1.lower() == closure[:len(kmer1)].lower(): 
                        print("kmers mismatch:",kmer1,closure,len(kmer1),len(closure), file=sys.stderr)
                        raise Exception 
                    outfasta.write( closure[len(kmer1):] )
                    nibble_off = len(kmer2)
                    last_patch = closure
                    nibbled_seq = closure[-len(kmer2):]
                elif (end,last_end) in gapclose_data:
                    kmer2,kmer1,closure= gapclose_data[end,last_end]
                    kmer1 = rc(kmer1)
                    kmer2 = rc(kmer2)
                    closure = rc(closure)

                    last_patch = closure
                    if not kmer1 == closure[:len(kmer1)]: 
                        print("kmers mismatch:",kmer1,closure,len(kmer1),len(closure), file=sys.stderr)
                        raise Exception 

                    outfasta.write( closure[len(kmer1):] )
                    nibble_off = len(kmer2)
                    nibbled_seq = closure[-len(kmer2):]

                    pass
                else:
                    outfasta.write( "N" * int(x-last_x)  )
                    nibbled_seq = ""
                    last_patch = ""

#            print(seqs[contig][:500], file=sys.stderr)
            if p=="5":
                outfasta.write(seqs[contig][nibble_off:])
                trailing_kmer = seqs[contig][-kmer:]
                if not seqs[contig][:nibble_off].lower()==nibbled_seq.lower():
                    print("inconsistent gap closure:",contig,nibble_off,seqs[contig][:nibble_off].lower(),"neq",nibbled_seq.lower(), file=sys.stderr)
                    raise Exception 
                    #exit(1)
            else:
                outfasta.write(rc(seqs[contig])[nibble_off:])
                trailing_kmer = rc(seqs[contig])[-kmer:]
                if not rc(seqs[contig])[:nibble_off].lower()==nibbled_seq.lower():

                    print("inconsistent gap closure:",contig,nibble_off,rc(seqs[contig])[:nibble_off],"neq",nibbled_seq,"rc", file=sys.stderr)
                    raise Exception 
                    #exit(1)

        else:
            pass
#            print("#lc:",contig,p,int(x-last_x),len(seqs[contig]), file=sys.stderr)
    #        print "N" * int(x-last_x)
        last_c=contig
        last_x=x
        last_s=scaffold
        last_end=end

    outfasta.flush()
    outfasta.close()

#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

import sys
import bisect


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

class MapperLite:
    def __init__(self):
        self.ocontig_contigs ={}
        self.contig_ocontig  ={}
        self.contig_scaffold ={}
        self.contigx         ={}
        self.contig_strand   ={}
        self.scaffolds       ={}
        self.contig_length   ={}
        self.scaffold_gaps   ={}
        self.gap2ends        ={}
        self.scaffold_contigs={}
        self.scaffold_length ={}

    def init_from_layout_list(self,layout_lines,scaffold_lengths):
        ocontig_contigs = self.ocontig_contigs
        contig_scaffold = self.contig_scaffold
        contigx         = self.contigx
        contig_strand   = self.contig_strand
        scaffolds       = self.scaffolds
        contig_length   = self.contig_length
        scaffold_gaps   = self.scaffold_gaps
        gap2ends        = self.gap2ends
        scaffold_contigs = self.scaffold_contigs
        self.scaffold_length = dict(scaffold_lengths)

        for i in range(len(layout_lines)-1):
            scaffold1,ocontig1,base1,end1,x1,l21 = self.layout_lines[i  ]
            scaffold2,ocontig2,base2,end2,x2,l22 = self.layout_lines[i+1]
            scaffolds[scaffold1]=1

            if not scaffold1==scaffold2: continue
            if ocontig1==ocontig2 and base1==base2:
                contig = ocontig1+"_"+str(base)
                self.contig_ocontig[contig]=ocontig1
                scaffold_contigs[scaffold] = scaffold_contigs.get(scaffold,[])+[contig]
                contig_length[contig]=l21
                contig_scaffold[contig]=scaffold1
                ocontig_contigs[ocontig1] = ocontig_contigs.get(ocontig,[])+[(int(base1)-1,int(base1)-1+contig_length,contig)]
                if end1=="5": 
                    contig_strand[contig]=1
                    contigx[contig]=x1
                else:
                    contig_strand[contig]=-1
                    contigx[contig]=x2


            else:
                scaffold_gaps[scaffold1]=scaffold_gaps.get(scaffold1,[]) + [(x1,x2)]
                gap2ends[(scaffold,x1,x2)] = ("{}_{}.{}".format(ocontig1,base1,end1),"{}_{}.{}".format(ocontig2,base2,end2))

        for ocontig in ocontig_contigs.keys():
            ocontig_contigs[ocontig].sort()
            

    def load_layout(self,fileh):
        ocontig_contigs = self.ocontig_contigs
        contig_scaffold = self.contig_scaffold
        contigx         = self.contigx
        contig_strand   = self.contig_strand
        scaffolds       = self.scaffolds
        contig_length   = self.contig_length
        scaffold_gaps   = self.scaffold_gaps
        gap2ends        = self.gap2ends
        scaffold_contigs = self.scaffold_contigs
        scaffold_length  = self.scaffold_length
        if 'readline' in dir(fileh):
            fh=fileh
        else:
            fh=open(fileh,"rt")
        
        last_contig = -1
        last_scaffold= -1

        while True:
            l = fh.readline()
            if not l: break

            if l[0]=="#": continue
            c=l.strip().split()
            if not c[0]=="p:": continue

            scaffold,end = c[1],c[2]
            scaffold_length[scaffold]=int(c[6])
            this_contig_length=int(c[7])
            x=int(float(c[3]))
            scaffolds[scaffold]=1
            contig = end[:-2]
            npp=contig.split("_")
            ocontig,base="_".join(npp[:-1]),npp[-1]
            self.contig_ocontig[contig]=ocontig
            base=int(base)

            if not last_contig==contig:
                 contig_strand[contig]=1 if (end[-1]=="5") else -1

            if end[-1]=="5": 
                contigx[contig]=x
                ocontig_contigs[ocontig] = ocontig_contigs.get(ocontig,[])+[(base-1,base-1+this_contig_length,contig)]
                scaffold_contigs[scaffold] = scaffold_contigs.get(scaffold,[])+[contig]

                ocontig_contigs[ocontig].sort()

            length=int(c[7])
            contig_length[contig]=length
            blast_coord = int(c[5])

            if last_scaffold==scaffold and not last_contig==contig: 
                scaffold_gaps[scaffold]=scaffold_gaps.get(scaffold,[]) + [(lastx,x)]
                gap2ends[(scaffold,lastx,x)] = (last_end,end)
            elif not last_scaffold==scaffold:
                pass
    #            print "gaps:",scaffold,scaffold_gaps.get(last_scaffold,[])
            blast_chr = c[4]

            contig_scaffold[contig] =scaffold
            contig_scaffold[ocontig]=scaffold

            bh=eval(" ".join(c[8:]))
            if bh:
                unc =abs(length- abs(int(bh[3])-int(bh[4])))

                if bh[2]=="+" and end[-1]=="5": blast_coord -= int(bh[5])-1
                if bh[2]=="+" and end[-1]=="3": blast_coord += length-int(bh[6])-1
                if bh[2]=="-" and end[-1]=="5": blast_coord += int(bh[5])-1
                if bh[2]=="-" and end[-1]=="3": blast_coord -= length-int(bh[6])-1

            else:
                unc=0

            last_scaffold=scaffold
            last_contig=contig
            lastx=x
            last_end=end

    def mapCoord(self,ocontig,x,y):        
        i=0
        for b,limit,contig in self.ocontig_contigs.get(ocontig,[])[::-1]:
            if b<=x and x<limit:
                if self.contig_strand[contig]==1:
#                    if self.contigx[contig]+(x-b) < 0 :
#                        print("coordinate mapping weirdness 1:",x,self.contigx[contig]+(x-b),contig,self.contigx[contig],b,i,len(self.ocontig_contigs.get(ocontig,[])))
#                        raise Exception 
#                    if self.contigx[contig]+(y-b) < 0 :
#                        print("coordinate mapping weirdness 2:",y,self.contigx[contig]+(y-b),contig,self.contigx[contig],b,i,len(self.ocontig_contigs.get(ocontig,[])))
#                        raise Exception 
                    return self.contig_scaffold[contig],self.contigx[contig]+(x-b),self.contigx[contig]+(y-b),self.contig_strand[contig],contig
                else:
#                    if self.contigx[contig]-(x-b) < 0 :
#                        print("coordinate mapping weirdness 3:",x,self.contigx[contig]-(x-b),contig,self.contigx[contig],b,i,len(self.ocontig_contigs.get(ocontig,[])))
#                        raise Exception 
#                    if self.contigx[contig]-(y-b) < 0 :
#                        print("coordinate mapping weirdness 4:",y,self.contigx[contig]-(y-b),contig,self.contigx[contig],b,i,len(self.ocontig_contigs.get(ocontig,[])))
#                        raise Exception 
                    return self.contig_scaffold[contig],self.contigx[contig]-(x-b),self.contigx[contig]-(y-b),self.contig_strand[contig],contig
            i+=1
        return [False] * 5

    def hits_gaps(self,scaffold,x,y,debug=False):
        gaps=[]
        a=self.scaffold_gaps.get(scaffold,[])
        i = bisect.bisect_left(a,(x,y))
#        if i>=0 and i<len(a):
#            return len(a),i,a[i],(x,y)
        j=i
        if j>=len(a): return([])
        while j>0 and a[j][1]>=x: j-=1

        while j<len(a) and a[j][0]<y:
            if debug: print("# hits overlap:",j,a[j],(x,y),y-x,pairs_overlap((x,y),a[j]))
            if pairs_overlap((x,y),a[j]):
                gaps.append(a[j])
            j=j+1
        return gaps


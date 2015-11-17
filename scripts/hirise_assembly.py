#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

import os
import re
import pysam
from bamtags import BamTags
import gc
import random
import bisect
import numpy as np

def parse_pline(l,scaffold_lengths):
    c=l.strip().split()
    m=re.match("(.*)_(\d+)\.([35])",c[2])
    scaffold_lengths[c[1]]=int(c[6])
    return( c[1],m.group(1),int(m.group(2)),m.group(3),int(c[3]),int(c[7]) )

def parse_pline2(l,scaffold_lengths):
    c=l.strip().split()
    base=int(c[3])
    x=int(c[5])
    scaffold_lengths[c[1]]=max(scaffold_lengths.get(c[1],0),x)
    return( c[1],c[2],base,c[4],x,int(c[6]) )

def aln2bins(aln,binsize):        
    return( ( aln.tid,int(aln.pos/binsize) ),( aln.rnext,int(aln.pnext/binsize)))


class GrowingNPArray:
    def __init__(self,dtype):
        self.data = np.zeros((10000,),dtype=dtype)
        self.N=10000
        self.n=0
        self.dtype = dtype

    def append(self,d):
        if self.n==self.N:
            self.N *=2
            dd = np.zeros((self.N,),dtype=self.dtype)
            dd[:self.n] = self.data
            self.data = dd
        self.data[self.n]=d
        self.n+=1

    def finalize(self):
        dd = np.zeros((self.n,),dtype=self.dtype)
        dd[:self.n] = self.data[:self.n]
        return dd


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
                    
class HiriseAssembly:
    """This class is designed to encapsulate a Hirise Assembly."""

    def __init__(self,options={}):
        """Options currently used are:  bams, bam, datamodel, layout"""
        self.bams=[]
        if options.get("bams"): self.bams+=list(options["bams"])
        if options.get("bam") : self.bams +=[options["bam"]]
        self.model_params=False

        self.binsize=10000
        self.bam_objects=False
        self.shotgun_bam_objects=False

        if "bam" in options or "bams" in options:
            self.load_ocontig_info()
        
        if "datamodel" in options and options["datamodel"]: 
            self.set_modelparams(options["datamodel"])

        self.layout_lines=[]
        if options.get("layout"):
            self.load_playout(options["layout"])
        self.masked_regions={}
        self.check_bams()
        self.ocontig_fa=False
        self.ocontig_fasta=False
        if "shotgun" in options: 
            self.shotgun = options["shotgun"]
        else:
            self.shotgun = []
        self.contigx=False

    def get_seq(self,contig):
        if not self.ocontig_fasta:
            raise Exception 
        
        if not os.path.exists(self.ocontig_fasta+".fai"):
            pysam.faidx(self.ocontig_fasta)
        if not self.ocontig_fa:
            self.ocontig_fa = pysam.Fastafile(self.ocontig_fasta)

        m=re.match("(.*)_(\d+)",contig)
        ocontig=m.group(1)
        x=int(m.group(2))
        y=x+self.contig_length(contig)
        leng=self.ocontig_lengths[ocontig]
        #return self.ocontig_lengths[m.group(1)]


        #dummy,scaffold,x,y,leng = l.strip().split()
        #x=int(x)
        #y=int(y)
        #leng=int(leng)
#        y=min(y,leng)
        seq = self.ocontig_fa.fetch(ocontig,x-1,y).decode()
        return(seq)

    def contig_length(self,contig):
        if contig in self.contig_lengths:
            return self.contig_lengths[contig]
        else:
            m=re.match("(.*)_(\d+)",contig)
            return self.ocontig_lengths[m.group(1)]
        return 0

    def contig_end_to_region(self,contig_end_string,w):
        ((this_contig,this_base,this_end),) = re.findall("^(.*)_(\d+)\.([53])$",contig_end_string)
        print(this_contig,this_base,this_end)
        this_base = int(this_base)
#        this_end = int(this_end)
        for scaffold,contig,base,end,x,l in self.layout_lines:
            slen = self.scaffold_lengths[scaffold]
            if this_contig == contig and this_base == base and this_end == end:
                return(scaffold,max(0,x-w),min(slen,x+w),(x<slen-x))
        print("so such end found",contig_end_string)
        raise Exception

    def index_ocontigs(self):
        self.ocontig_bases  = {}
        self.contig_lengths = {}
        self.contig_scaffold= {}
        self.load_ocontig_info()
#        self.scaffold_ocontigs={}
        
        ocb=self.ocontig_bases
        for scaffold,contig,base,end,x,l in self.layout_lines:

#            if not scaffold in self.scaffold_ocontigs:
#                self.scaffold_ocontigs[scaffold]={}
#            self.scaffold_ocontigs[scaffold][self.contig_ocontig[contig]]=1

            self.contig_scaffold[contig,base]=scaffold
            if end in [5,"5"] and not base==1:
                ocb[contig] = ocb.get(contig,[])+[ base ]
            if end in [5,"5"]:  
                self.contig_lengths[ contig+"_"+str(base) ]=l
                self.contig_lengths[ contig,base ]=l
        for b in ocb.keys():
            ocb[b].sort()
            #print("ocb:",b,ocb[b])
        self.total_ocontig_length = sum(self.ocontig_lengths.values())
        self.rc_buffer=[]
        x=0
        for i in self.ocontigs_iter():
            self.rc_buffer.append((x,i))
            x+=self.ocontig_lengths[i]

    def broken_coord(self,ocontig,x,tuples=False):
        if tuples:
            if not ocontig in self.ocontig_bases:
                if (not (ocontig,1) in self.contig_lengths ) or x < self.contig_lengths[(ocontig,1)]: 
#                    print(ocontig,"1")
                    return((ocontig,1), x)
                else:
#                    print(ocontig,"2")
                    return((ocontig, self.contig_lengths[(ocontig,1)]), x)
            for b in self.ocontig_bases[ocontig][::-1]:
                if b<x:
                    new_contig = (ocontig,b)
                    #print("zxzx",ocontig,x,b,new_contig,self.contig_lengths.get(new_contig,False))
                    if (not new_contig in self.contig_lengths) or (x-b)<self.contig_lengths[new_contig]:
#                        print(ocontig,"3")
                        return((ocontig,b), x-b)
                    else:
                        return (False,False)
            new_contig = (ocontig,1)
            if (not new_contig in self.contig_lengths) or x < self.contig_lengths[new_contig]:
#                print(ocontig,"4",new_contig,new_contig in self.contig_lengths)
                return(new_contig, x)
            else:
                return (False,False)
            

        else:
            if not ocontig in self.ocontig_bases:
                if (not ocontig+"_1" in self.contig_lengths ) or x < self.contig_lengths[ocontig+"_1"]: 
                    return(ocontig+"_1", x)
                else:
                    return(ocontig+"_"+str( self.contig_lengths[ocontig+"_1"] ), x)
            for b in self.ocontig_bases[ocontig][::-1]:
                if b<x:
                    new_contig = ocontig+"_"+str(b)
                    #print("zxzx",ocontig,x,b,new_contig,self.contig_lengths.get(new_contig,False))
                    if (not new_contig in self.contig_lengths) or (x-b)<self.contig_lengths[new_contig]:
                        return(ocontig+"_"+str(b), x-b)
                    else:
                        return (False,False)
            new_contig = ocontig+"_1"
            if (not new_contig in self.contig_lengths) or x < self.contig_lengths[new_contig]:
                return(new_contig, x)
            else:
                return (False,False)

    def merge_masked_regions(self,debug=False):
        mr = self.masked_regions
        for ocontig in mr.keys():
            new_regions=[]
            edges = []
            for x,y in mr[ocontig]:
                edges.append((x,1))
                edges.append((y,-1))
            edges.sort(key=lambda x:(x[0],-x[1]))
            rs=0
            state=0
            for x,b in edges:
                rs+=b
                if state==0 and rs>0:
                    start=x
                    state=1
                elif state==1 and rs==0:
                    end=x
                    new_regions.append((start,end))
                    state=0
            if debug and not (len(mr[ocontig])==len(new_regions)):
                print("#mmr",ocontig,len(mr[ocontig]),len(new_regions),mr[ocontig],new_regions)
            mr[ocontig]=new_regions

    def add_mask_regions(self,segments=[],filename=False):
        masked_regions = self.masked_regions
        if filename:
            for line in open(filename):
                c=line.strip().split()
                ocontig,x,y = c[0],int(c[1]),int(c[2])
                if not ocontig in masked_regions: masked_regions[ocontig]=[]
                masked_regions[ocontig].append((x,y))
        for ocontig,x,y in segments:
            if not ocontig in masked_regions: masked_regions[ocontig]=[]
            masked_regions[ocontig].append((x,y))
            
        for oc in masked_regions.keys():
            masked_regions[oc].sort()

    def setup_mapper(self,debug=False):
#        debug=True
        self.contigx={}
        self.contigst={}
        cex={}
        for scaffold,contig,base,end,x,l in self.layout_lines:
            if end in [5,"5"]:
                self.contigx[contig,base] = x
            cex[contig,base,end]=x
        for contig,base in self.contigx.keys():
            if cex[contig,base,"5"] < cex[contig,base,"3"]: 
                self.contigst[contig,base]=1
            else: 
                self.contigst[contig,base]=-1
#        if debug:
#            for contig,base in self.contigx.keys():
#                print("m:",contig,base,self.contigx[contig,base],self.contigst[contig,base])


    def scaffolds_iter(self,debug=False):
        if not self.layout_lines or len(self.layout_lines)==0:
            self.make_trivial_layout_lines(debug=debug)
            self.validate()

        last_scaffold=False
        for scaffold,contig,base,end,x,l in self.layout_lines:
            if (not last_scaffold) or not scaffold==last_scaffold:
                last_scaffold=scaffold
#                print("#",scaffold)
                yield(scaffold)

    def scaffold_coord(self,c,x,debug=False):
        if not self.contigx: self.setup_mapper(debug=debug)
        y=self.contigx[c]+self.contigst[c]*x
        if debug: print("sc:",c,x,y)
        return y

    def chicago_support_scan(self,scaffold,minsupport=0.0,mapq=10,debug=False,minx=0,maxx=1e6,gpf=False,logfile=False):

        import chicago_edge_scores as ces
        ces.set_exp_insert_size_dist_fit_params(self.model_params)
        model=ces.model

        links={}
        self.get_scaffold_links(scaffold,skipI=False,mapq=mapq,links=links,debug=debug,tuples=True)
#        model = ces.
        pairs=[]
        edges=[]
        buff=100
        nb=[]
        joins=[]
        t1=50.0
        t2=minsupport
        slen=self.ocontig_lengths[scaffold]
        for c1,c2 in links.keys():
            if debug: print(c1,c2,links[c1,c2])
            for x,y in links[c1,c2]:
                pairs.append( ( self.scaffold_coord(c1,x,debug=debug),self.scaffold_coord(c2,y,debug=debug)) )
        for a,b in pairs:               
            a,b=min(a,b),max(a,b)
            if debug:
                print("z",a,b)


            if (b-a)>minx and (b-a)<maxx:
                if gpf: gpf.write("{}\t{}\t{}\t{}\n".format(0.5*(a+b),0.5*(b-a),w,z))
                ll=model.lnF(b-a)

                edges.append( tuple([a+buff         , ll,c1,"a"]) )
                edges.append( tuple([b-buff         ,-ll,c2,"b"]) )

                nb.append( tuple([a+buff         , 1]) )
                nb.append( tuple([b-buff         ,-1]) )

        if gpf: gpf.write("\n\n\n")
        edges.sort()
        nb.sort()
        rs=0
        n=0
        tripped=False
        last=0
        state=0
        stretch_start=0
        low_point=0.0
        ji=0
        gap_scores={}
        gap_coverages={}
        max_t1=0
        min_t1=slen
        breakers=[]
        
        for i in range(len(edges)):
            rs+=edges[i][1]
            n+=nb[i][1]
            x=edges[i][0]

            try:
                score=model.cutScore(slen,x,n,rs,rangeCutoff=maxx,minSep=minx)
            except Exception as e:
                print(edges)
                print(i,edges[i],nb[i],rs,n,x,scaffold,slen,len(edges))
                raise e


            if score>t1:
                tripped=True
                last_t1=x
                min_t1 = min(min_t1,x)
                max_t1 = max(max_t1,x)
            if tripped and score<t2 and state==0:
                stretch_start=x
                state=1
                low_point =score
                minx = x
            if state==1 and score>t2:
                breakers.append((scaffold,stretch_start,x,minx,low_point,slen))
#                if logfile: logfile.write("{} {} {} {} {} {}\n".format(scaffold,stretch_start,x,minx,low_point,slen))
                state=0
            if state==1:
                if score<low_point:
                    low_point=score
                    minx=x
            if debug: print("dd:",scaffold,edges[i][0],score,rs,state,x,stretch_start,score,n,edges[i][2])
            if gpf: gpf.write("{}\t{}\t{}\t{}\n".format(edges[i][0],score,rs,n))

            last=edges[i][0]
            
        for scaffold,stretch_start,x,minx,low_point,slen in breakers:
            if stretch_start < last_t1:
#                if logfile: logfile.write("{} {} {} {} {} {} {} {}\n".format(scaffold,stretch_start,x,minx,low_point,slen,min_t1,max_t1))
                if logfile: logfile.write("{} {} {} {} {} {}\n".format(scaffold,stretch_start,x,minx,low_point,slen))

    def make_trivial_layout_lines(self,debug=False):
        self.layout_lines=[]
        #print(self.ocontig_lengths)
        if len(self.ocontig_lengths)==0:
            print("#load_contig_info")
            self.load_ocontig_info()
        i=1
        for ocontig in self.ocontigs_iter():
            scaffold=ocontig
#            if debug: print(scaffold,ocontig,1,"5",0,self.ocontig_lengths[ocontig],sep="\t")
            self.layout_lines.append((scaffold,ocontig,1,"5",0,self.ocontig_lengths[ocontig])) #
            self.layout_lines.append((scaffold,ocontig,1,"3",self.ocontig_lengths[ocontig],self.ocontig_lengths[ocontig])) #
            self.scaffold_lengths[scaffold] = self.ocontig_lengths[ocontig]
            self.scaffold_ocontigs[scaffold] = self.scaffold_ocontigs.get(scaffold,[])+[ocontig]
            i+=1
        self.index_ocontigs()
        self.setup_mapper()

    def get_scaffold_links(self,Tscaffold,skipI=False,mapq=10,links={},debug=False,tuples=False):
        if not self.layout_lines or len(self.layout_lines)==0:
            self.make_trivial_layout_lines(debug=debug)
            self.validate()
        contigs={}
        ocontigs={}
        if not tuples:
            for scaffold,contig,base,end,x,l in self.layout_lines:
                if scaffold==Tscaffold and end in [5,"5"]:
                    contigs[ contig+"_"+str(base) ]=1  
                    ocontigs[contig]=1
        else:
            for scaffold,contig,base,end,x,l in self.layout_lines:
                if scaffold==Tscaffold and end in [5,"5"]:
                    contigs[ contig,base ]=1  
                    ocontigs[contig]=1
                
        if debug: print(ocontigs,contigs)
        self.get_links(list(ocontigs.keys()),skipI,mapq,links,contigs,tuples=tuples)

    def segments_iter(self):
        for l in self.layout_lines:
            scaffold1,contig1,base1,end1,x1,l21 =l
            if end1 in [5 , "5"]:
                #cc=self.contig_length(contig1+"_"+str(base1))
                yield( contig1,base1,base1+l21,self.ocontig_lengths[contig1] )

    def contigs_iter(self):
        for l in self.layout_lines:
            scaffold1,contig1,base1,end1,x1,l21 =l
            if end1 in [5 , "5"]:
                yield( contig1 +"_"+str(base1) )

    def ocontigs_iter(self):
        for ocontig in self.ocontig_lengths.keys():
            yield(ocontig)

    def load_ocontig_info(self):
        import pysam

        self.ocontig_lengths={}
        if not self.bam_objects:
            self.bam_objects={}
            for b in self.bams:
                self.bam_objects[b] =  pysam.Samfile(b,"rb") 
        for bam in self.bam_objects.values():
            h=bam.header
            seqs=h['SQ']
            for s in seqs:
                if s['SN'] in self.ocontig_lengths:
                    if not  self.ocontig_lengths[s['SN']]==s['LN'] :
                        raise Exception 
                else:
                    self.ocontig_lengths[s['SN']]=s['LN']

    def random_window(self,wlen=1000):
        r=int(random.random()*self.total_ocontig_length)
        x=bisect.bisect_right(self.rc_buffer,(r,"a"))
        while r+x > self.total_ocontig_length or (x<len(self.rc_buffer)-1 and r+wlen>self.rc_buffer[x+1][0]):
            r=int(random.random()*self.total_ocontig_length)
            x=bisect.bisect_right(self.rc_buffer,(r,"a"))
#            print( r,self.rc_buffer[x],self.rc_buffer[x+1],r-self.rc_buffer[x][0])
        x-=1
#        print( r,self.rc_buffer[x],self.rc_buffer[x+1],r-self.rc_buffer[x][0])
        return( self.rc_buffer[x][1], r-self.rc_buffer[x][0], r-self.rc_buffer[x][0]+wlen )

    

    def top_other_contigs(self,reference,mapq=10,n=2,debug=False,bins=False,binx=0):
        nother={}
        if not self.bam_objects:
            self.bam_objects={}
            for b in self.bams:
                self.bam_objects[b] = pysam.Samfile(b,"rb")
        if debug: print(self.bam_objects)

        if bins:
            region="{}:{}-{}".format(reference,binx*self.binsize,(binx+1)*self.binsize)
            for bam in self.bam_objects.values():
                for aln in bam.fetch(region=region,until_eof=True):
                    if aln.is_duplicate :              continue
                    if aln.mapq < mapq               : continue
                    if BamTags.mate_mapq(aln) < mapq : continue
                    my_bin,rbin = aln2bins(aln,self.binsize)
                    if not my_bin == rbin : 
                        nother[rbin] = nother.get(rbin,0)+1
        else:
            for bam in self.bam_objects.values():
                for aln in bam.fetch(reference=reference,until_eof=True):
                    if aln.is_duplicate :              continue
                    if aln.mapq < mapq               : continue
                    if BamTags.mate_mapq(aln) < mapq : continue
                    if not aln.rnext == aln.tid : 
                        nother[aln.rnext] = nother.get(aln.rnext,0)+1
 
        kk=list(nother.keys())
        kk.sort(key=lambda x: nother[x],reverse=True)
        if debug: print(kk[:2])
        return(tuple(kk[:2]))

    def window_stats(self,reference,xa,xb,mapq=10,debug=False,bins=False):
        if not self.bam_objects:
            self.bam_objects=[]
            for b in self.bams:
                self.bam_objects.append( pysam.Samfile(b,"rb") )

        if not self.shotgun_bam_objects:
            self.shotgun_bam_objects=[]
            for b in self.shotgun:
                self.shotgun_bam_objects.append( pysam.Samfile(b,"rb") )
        
        binx = 0
        if bins:
            binx = int(xa/self.binsize)

        top2 = self.top_other_contigs(reference,mapq=mapq,debug=debug,bins=bins,binx=binx)
        region = "{}:{}-{}".format(reference,max(1,xa-1000),xb+1000)
        n_aligned=0
        n_hiq=0
        templates={}
        n_lowq   =0
        occ ={}
        occ2={}
        for bam in self.bam_objects.values():
            if debug: print("#",region,xa,xb,bam)
            for aln in bam.fetch(region=region):
                if aln.pos > xb:                   continue
                if aln.pos < xa:                   continue
                if aln.is_duplicate :              continue
                n_aligned+=1
                if aln.mapq < mapq               : continue
                if BamTags.mate_mapq(aln) < mapq : continue
                templates[aln.query_name]=templates.get(aln.query_name,0)+1
                if bins:
                    my_bin,rbin = aln2bins(aln,self.binsize)
                    if not my_bin==rbin:
                        occ[rbin] = occ.get(rbin,0)+1
                        if not rbin in top2:
                            occ2[rbin] = occ2.get(rbin,0)+1
                            
                else:
                    if (not aln.rnext == aln.tid):
                        occ[aln.rnext] = occ.get(aln.rnext,0)+1
                        if (not aln.rnext in top2) : 
                            occ2[aln.rnext] = occ2.get(aln.rnext,0)+1
                n_hiq+=1
            if debug: print("#",n_hiq,n_aligned)
        n_both   = len([ i for i in templates.keys() if templates[i]>1  ])
        n_other  = len(occ.keys())
        n_other2 = len([ i for i in occ.keys() if occ[i]>1  ])
        n_other3 = len([ i for i in occ.keys() if occ[i]>2  ])

        n_otherB  = len(occ2.keys())
        n_otherB2 = len([ i for i in occ2.keys() if occ[i]>1  ])
        n_otherB3 = len([ i for i in occ2.keys() if occ[i]>2  ])

        n_shotgun=0
        for bam in self.shotgun_bam_objects:
            for aln in bam.fetch(region=region):
                if aln.pos > xb: continue
                if aln.is_duplicate :              continue
                if aln.mapq < mapq               : continue
                n_shotgun+=1
        
        print("\t".join(map(str,[region,n_aligned,n_hiq,n_both,n_other,n_other2,n_other3,n_otherB,n_otherB2,n_otherB3,n_shotgun])))

    def chicago_promiscuity_mask(self,ocontig,mapq=10,w=1000,minlinks=2,maxothers=3,outfile=False,bins=False,debug=False):
        if not self.bam_objects:
            self.bam_objects=[]
            for b in self.bams:
                self.bam_objects.append( pysam.Samfile(b,"rb") )
#        dtype=[('x',int),('t',int),('y',int),('b',int),('r',int)]
        
        dtype=[('x',int),('t',int),('y',int)]
        buffer = GrowingNPArray(dtype=dtype)

        i=0
        tid=0
        for bam in self.bam_objects:                
            if debug: print("#",ocontig,bam)
            tid = bam.gettid(ocontig)
            for aln in bam.fetch(reference=ocontig):
                if aln.rnext<0: continue
                if aln.mapq < mapq               : continue
                if BamTags.mate_mapq(aln) < mapq : continue
                if aln.is_duplicate : continue
#                buffer.append((aln.pos,aln.rnext,aln.pnext,i,aln.tid))
                buffer.append((aln.pos,aln.rnext,aln.pnext))
            i+=1

        buffer = buffer.finalize()
        buffer = np.sort(buffer,order='x',kind="mergesort")

        bad_edges=[]


        def get_contig(x,b):
            return x[1]

        def get_contig_bin(x,b):
            return(x[1]+(int(x[2]/b)/1.0e6))

        if bins:
            get_bin = get_contig_bin
        else:
            get_bin = get_contig

        bc={}
        i=0
        j=0
        N=len(buffer)
        while i < N:
            mybin = get_bin((0,tid,buffer[i][0]),self.binsize)
            while j<N and buffer[j][0]<buffer[i][0]+w: 
                jbin = get_bin(buffer[j],self.binsize)
                bc[ jbin ] = bc.get(jbin,0)+1
                j+=1
                
            n_in_window = j-i
            
            nw=len([ ii for ii in bc.keys() if bc[ii]>=minlinks and not ii==mybin ])  #n_other_windows(buffer[i:j],minlinks,bins,self.binsize)

#            if debug: print(tid,buffer[i][0],n_in_window,nw,mybin,[ (ii,bc[ii]) for ii in bc.keys() if bc[ii]>=minlinks and not ii==mybin ],[ get_bin(x,self.binsize) for x in buffer[i:j] ])
            if debug: print(tid,buffer[i][0],n_in_window,nw,mybin,[ (ii,bc[ii]) for ii in bc.keys() if bc[ii]>=minlinks and not ii==mybin ])


            if nw > maxothers:
                bad_edges.append((buffer[i  ][0], 1))
                bad_edges.append((buffer[j-1][0],-1))

            ibin = get_bin(buffer[i],self.binsize)
            bc[ ibin ] = bc[ibin]-1
            if bc[ibin]==0: bc.pop(ibin)
            i+=1

        bad_edges.sort()
        if debug:
            for be in bad_edges:
                print(e)

        rs=0
        state=0
        for x,dy in bad_edges:
            rs+=dy
            if state==0 and rs>0:
                startx=x
                state=1
            elif state==1 and rs==0:
                outfile.write("{} {} {}\n".format(ocontig,startx,x))
                if debug: print("#{} {} {}\n".format(ocontig,startx,x))
                state=0
        outfile.flush()

        
    def chicago_read_density_mask(self,ocontigs,mapq=10,outfile=False,w=1000,cutoff=100,debug=False,shotgunT=False):

        if not self.bam_objects:
            self.bam_objects=[]
            for b in self.bams:
                self.bam_objects.append( pysam.Samfile(b,"rb") )

        if not self.shotgun_bam_objects:
            self.shotgun_bam_objects=[]
            for b in self.shotgun:
                self.shotgun_bam_objects.append( pysam.Samfile(b,"rb") )

        for ocontig in ocontigs:
            locations=[]
            shotgun_locations=[]
#            gc.disable()
            for bam in self.bam_objects:                
                if debug: print("#",ocontig,bam)
                for aln in bam.fetch(reference=ocontig):
#                    if debug: print("#",aln.pos,aln.is_duplicate,aln.mapq,BamTags.mate_mapq(aln),aln.query_name,len(locations))
                    if aln.rnext<0: continue
                    if aln.mapq < mapq               : continue
                    if BamTags.mate_mapq(aln) < mapq : continue
                    if aln.is_duplicate : continue
                    locations.append(aln.pos)
                if debug: print("#",len(locations))
            if shotgunT:
                for bam in self.shotgun_bam_objects:                
                    if debug: print("#",ocontig,bam)
                    for aln in bam.fetch(reference=ocontig):
                        if aln.mapq < mapq  : continue
                        if aln.is_duplicate : continue
                        shotgun_locations.append(aln.pos)
                    if debug: print("#",len(locations),len(shotgun_locations))

#            gc.enable()
            locations.sort()
            shotgun_locations.sort()

            i=0
            j=0
            bad_edges=[]

            while i<len(locations):
                while j<len(locations) and locations[j]<locations[i]+w: j+=1
                n_in_window = j-i
                if debug:
                    print("chicago_density",ocontig,locations[i],locations[j-1],j-i,locations[j-1]-locations[i],locations[i:j])
                if n_in_window > cutoff:
    #                outfile.write("{} {} {}\n".format(ocontig,locations[i],locations[j-1]))
                    bad_edges.append((locations[i],1)  )
                    bad_edges.append((locations[j-1],-1)  )
                i+=1

            if shotgunT:
                i=0
                j=0
                while i<len(shotgun_locations):
                    while j<len(shotgun_locations) and shotgun_locations[j]<shotgun_locations[i]+w: j+=1
                    n_in_window = j-i
                    if debug:
                        print("shotgun_density",ocontig,shotgun_locations[i],shotgun_locations[j-1],j-i,shotgun_locations[j-1]-shotgun_locations[i])
                    if n_in_window > shotgunT:
        #                outfile.write("{} {} {}\n".format(ocontig,locations[i],locations[j-1]))
                        bad_edges.append((shotgun_locations[i],1)  )
                        bad_edges.append((shotgun_locations[j-1],-1)  )
                    i+=1


            bad_edges.sort()

            rs=0
            state=0
            for x,dy in bad_edges:
                rs+=dy
                if state==0 and rs>0:
                    startx=x
                    state=1
                elif state==1 and rs==0:
                    outfile.write("{} {} {}\n".format(ocontig,startx,x))
                    state=0
        outfile.flush()


    def make_scaffold_mask(self,scaffold): #,mapper,mask):
        s_mask = {}
        segments={}
        slen = self.scaffold_lengths[scaffold]

        for ocontig in self.scaffold_ocontigs[scaffold]:
#            print(ocontig)
            for a,b in self.masked_regions.get(ocontig,[]):
                
                c1,x = self.broken_coord(ocontig,a,tuples=True)
                c2,y = self.broken_coord(ocontig,b,tuples=True)

                s1=self.contig_scaffold.get(c1)
                s2=self.contig_scaffold.get(c2)

                if c1==c2 and s1==s2 and s1==scaffold:

                    xx=self.scaffold_coord(c1,x)
                    yy=self.scaffold_coord(c2,y)
                    #print("#msm:",scaffold,ocontig,a,b,x,y,xx,yy,s1,s2,scaffold,sep="\t")

                    segments[ max(0,min(xx,yy)) , min( slen, max(xx,yy) )   ]=1
                elif (not c1==c2) and s1==s2 and s1==scaffold:
                    
                    xA = x
                    yA = self.contig_lengths[c1]
                    xx=self.scaffold_coord(c1,xA)
                    yy=self.scaffold_coord(c1,yA)
                    segments[ max(0,min(xx,yy)) , min( slen, max(xx,yy) )   ]=1

                    xB = 1
                    yB = y
                    xx=self.scaffold_coord(c2,xB)
                    yy=self.scaffold_coord(c2,yB)
                    segments[ max(0,min(xx,yy)) , min( slen, max(xx,yy) )   ]=1

                elif (not c1==c2) and s1==scaffold and not s2==scaffold:
                    
                    xA = x
                    yA = self.contig_lengths[c1]
                    xx=self.scaffold_coord(c1,xA)
                    yy=self.scaffold_coord(c1,yA)
                    segments[ max(0,min(xx,yy)) , min( slen, max(xx,yy) )   ]=1

                elif (not c1==c2) and s2==scaffold and not s1==scaffold:
                    
                    xB = 1
                    yB = y
                    xx=self.scaffold_coord(c2,xB)
                    yy=self.scaffold_coord(c2,yB)
                    segments[ max(0,min(xx,yy)) , min( slen, max(xx,yy) )   ]=1
 
        segments=list(segments.keys())
        segments.sort()
        #print("#msm:",scaffold,slen,segments,sep="\t")
        return(segments)

    def read_deserts(self,outfile,mapq=10,min_len=1000):
        if not self.bam_objects:
            self.bam_objects={}
            for b in self.bams:                
                self.bam_objects[b]= pysam.Samfile(b,"rb") 
        desert_edges=[]
        for bam in self.bam_objects.values():
            lastx=0
            lastr=-1
            for aln in bam.fetch(until_eof=True):
                if aln.mapq < mapq               : continue
                ref = aln.tid
                if not ref==lastr:
                    if lastr>-1:
                        lastname =bam.getrname(lastr)
                        lastl = self.ocontig_lengths[ lastname ]
                        if lastl - lastx > min_len:
                            desert_edges.append(( lastname,lastx  , 1 ))
                            desert_edges.append(( lastname,lastl  ,-1 ))
                    lastx=0
                deltax= aln.pos - lastx
                if deltax>min_len:
                    desert_edges.append(( bam.getrname(aln.tid),lastx  , 1 ))
                    desert_edges.append(( bam.getrname(aln.tid),aln.pos,-1 ))

                lastx=aln.pos
                lastr=aln.tid

        
        nbams = len( self.bams)
        state=0
        startx=0
        desert_edges.sort()
        y=0
        for contig,x,dy in desert_edges:

            y+=dy

            if state==0 and y == nbams:
                state=1
                startx=x
            elif state==1 and y<nbams:
                endx = x
                state=0
                if endx-startx>min_len:
                    print(contig,startx,endx,"readDesert",file=outfile)

            last_contig,last_x = contig,x


    def chicago_pairs(self,mapq=10,scaffold_contig_counts={},contig_scaffold_counts={},callback=False,bamfile=False,my_scaffolds=False):
        read_length = 100
        if not self.bam_objects:
            self.bam_objects={}
            for b in self.bams:                
                self.bam_objects[b]= pysam.Samfile(b,"rb") 
                


        if bamfile:
            bams_iter = iter( [ self.bam_objects[bamfile] ] )
        else:
            bams_iter = iter( self.bam_objects.values() )

        for bam in bams_iter:
#            print("#",bam)
            last_contig=False
            mri=0
            mr=[]
            last_scaffold=False
            for aln in bam.fetch(until_eof=True):
                if not aln.is_read1: continue
                track=False
                skip=False
                if aln.rnext<0:      continue
                if aln.mapq < mapq               : continue
                if BamTags.mate_mapq(aln) < mapq : continue
                if aln.is_duplicate : continue
                contig  = bam.getrname(aln.tid)
                contig2 = bam.getrname(aln.rnext)
                
                if not last_contig or not last_contig==contig:
                    if contig in self.masked_regions:
                        mr = self.masked_regions[contig]
                        mri=0
                    else:
                        mr=[]
                        mri=0
                c1,x = self.broken_coord(contig ,aln.pos  ,tuples=True)
                scaffold=self.contig_scaffold.get(c1)
                c2,y = self.broken_coord(contig2,aln.pnext,tuples=True)
                nscaffold=self.contig_scaffold.get(c2)

                if (not my_scaffolds==False) and (not ((scaffold in my_scaffolds) and (nscaffold in my_scaffolds))): continue

                if last_contig and not last_contig==contig:
#                    print("#done with contig",last_contig,list(contig_scaffold_counts.get(last_contig,{}).keys()))
                    for sc in contig_scaffold_counts.get(last_contig,{}).keys():  # for each hirise scaffold that this contig contributes a segment to:
                        if sc in scaffold_contig_counts:                   
                            if last_contig in scaffold_contig_counts[sc]:  # if we are still tracking the number of times we need to read this contig from a bam file for this scaffold
#                                print("##",sc,last_contig,scaffold_contig_counts[sc][last_contig])
                                scaffold_contig_counts[sc][last_contig] -= 1 # decrease the count by one
#                                print("##",sc,last_contig,scaffold_contig_counts[sc][last_contig])
                                if scaffold_contig_counts[sc][last_contig] ==0:
                                    del scaffold_contig_counts[sc][last_contig]
                            if len(scaffold_contig_counts[sc])==0:
#                                print("scaffold {} ready now",sc)
                                if callback: 
                                    callback(sc)
                                else:
                                    yield(sc,"done")
                                
                                del scaffold_contig_counts[sc]
#                                raise ScaffoldReady(sc)
#                            print("#s",last_contig,sc,scaffold_contig_counts.get(sc))
#                    print("#c",last_contig,contig_scaffold_counts.get(last_contig))
                last_contig = contig        
                last_scaffold = scaffold

                while mri<len(mr) and aln.pos > mr[mri][1]: mri+=1                        
                if mri<len(mr) and mr[mri][0]<aln.pos and aln.pos+read_length<mr[mri][1] : continue #this read maps to a masked region 

                if contig2 in self.masked_regions:
                    if [ (aa,bb) for aa,bb in self.masked_regions[contig2] if aa<aln.pnext and aln.pnext+read_length<bb ]: continue  #the pair maps to a masked region
                
                if not scaffold: continue
                if not nscaffold: continue

                xx=self.scaffold_coord(c1,x)
                yy=self.scaffold_coord(c2,y)

                #                yield(s1,self.scaffold_coord(c1,x),s2,self.scaffold_coord(c2,y))

                if self.scaffold_lengths[scaffold] < xx or  self.scaffold_lengths[nscaffold] < yy:
                    print("coordinate out of range error",xx,yy,scaffold,nscaffold,self.scaffold_lengths[scaffold], self.scaffold_lengths[nscaffold],contig,contig2,c1,c2,x,y)
                    raise Exception 


                yield( scaffold, nscaffold, xx, yy, c1, c2, aln.query_name )


    def chicago_pairs_for_scaffolds(self,mapq=10,callback=False,bamfile=False,scaffolds=[],contigs=False,minsep=0):
        read_length = 100
        if not self.bam_objects:
            self.bam_objects={}
            for b in self.bams:                
                self.bam_objects[b]= pysam.Samfile(b,"rb") 

        if bamfile:
            bams_iter = iter( [ self.bam_objects[bamfile] ] )
        else:
            bams_iter = iter( self.bam_objects.values() )


        ocontigs_to_hit = {}
        my_ocontigs=[]
        if contigs:
            my_ocontigs=contigs
        else:
            for scaffold in scaffolds:
                for oc in self.scaffold_ocontigs[scaffold]:
                    ocontigs_to_hit[oc]=1

            my_ocontigs = list(ocontigs_to_hit.keys())

        for bam in bams_iter:
#            print("#",bam)
            last_contig=False
            mri=0
            mr=[]
            last_scaffold=False
            for ocontig in my_ocontigs:
                for aln in bam.fetch(reference=ocontig):
                    if not aln.is_read1: continue
                    track=False
                    skip=False
                    if aln.rnext<0:      continue
                    if aln.mapq < mapq               : continue
                    if BamTags.mate_mapq(aln) < mapq : continue
                    if aln.is_duplicate : continue
                    contig  = bam.getrname(aln.tid)
                    contig2 = bam.getrname(aln.rnext)

                    if not last_contig or not last_contig==contig:
                        if contig in self.masked_regions:
                            mr = self.masked_regions[contig]
                            mri=0
                        else:
                            mr=[]
                            mri=0

                    if last_contig and not last_contig==contig:
                        pass
    #                    print("#done with contig",last_contig,list(contig_scaffold_counts.get(last_contig,{}).keys()))

                    c1,x = self.broken_coord(contig ,aln.pos  ,tuples=True)
                    scaffold=self.contig_scaffold.get(c1)
                    last_contig = contig        
                    last_scaffold = scaffold

                    while mri<len(mr) and aln.pos > mr[mri][1]: mri+=1                        
                    if mri<len(mr) and mr[mri][0]<aln.pos and aln.pos+read_length<mr[mri][1] : continue #this read maps to a masked region 

                    if contig2 in self.masked_regions:
                        if [ (aa,bb) for aa,bb in self.masked_regions[contig2] if aa<aln.pnext and aln.pnext+read_length<bb ]: continue  #the pair maps to a masked region

                    c2,y = self.broken_coord(contig2,aln.pnext,tuples=True)
                    nscaffold=self.contig_scaffold.get(c2)
                    if not scaffold: continue
                    if not nscaffold: continue

                    if not scaffold in scaffolds: continue
                    if not nscaffold in scaffolds: continue

                    xx=self.scaffold_coord(c1,x)
                    yy=self.scaffold_coord(c2,y)
                    if scaffold==nscaffold and abs(xx-yy)<minsep: continue

                    #                yield(s1,self.scaffold_coord(c1,x),s2,self.scaffold_coord(c2,y))

                    if self.scaffold_lengths[scaffold] < xx or  self.scaffold_lengths[nscaffold] < yy:
                        print("coordinate out of range error",xx,yy,scaffold,nscaffold,self.scaffold_lengths[scaffold], self.scaffold_lengths[nscaffold],contig,contig2,c1,c2,x,y)
                        raise Exception 

                    yield( scaffold, nscaffold, xx, yy, c1, c2, aln.query_name )

    def get_links(self,ocontigs,skipI=False,mapq=10,links={},contigs=False,raw=False,tuples=False,debug=False):
#        links={}

        if not self.bam_objects:
            self.bam_objects={}
            for b in self.bams:

                self.bam_objects[b]= pysam.Samfile(b,"rb") 

        for bam in self.bam_objects.values():
#            bam = pysam.Samfile(b,"rb") 

            for ocontig in ocontigs:
#                print("#",ocontig)
                mr=[]
                mri=0
                if ocontig in self.masked_regions: mr=self.masked_regions[ocontig]
                for aln in bam.fetch(reference=ocontig):
#                    if debug: print("#",aln.query_name,aln.pos,aln.pnext,aln.mapq,sep="\t")
                    if aln.rnext<0: continue
                    if aln.mapq < mapq               : continue
                    if BamTags.mate_mapq(aln) < mapq : continue
                    if aln.is_duplicate : continue

                    while mri<len(mr) and aln.pos > mr[mri][1]: mri+=1                        
                    if mri<len(mr) and mr[mri][0]<aln.pos and aln.pos<mr[mri][1] : 
                        if debug: print("#m1",ocontig,bam.getrname(aln.rnext),aln.pos,aln.pnext,mr[mri],aln.query_name,sep="\t")
                        continue #this read maps to a masked region 

                    if raw:
                        c1,x1 = ocontig,aln.pos
                    else:
                        c1,x1 = self.broken_coord(ocontig,aln.pos,tuples=tuples)
                        if not c1: continue

                    if contigs and not c1 in contigs: continue
                    ocontig2 = bam.getrname(aln.rnext)
                    
                    
                    if raw:
                        c2,x2 = ocontig2,aln.pnext
                    else:
                        c2,x2 = self.broken_coord(ocontig2,aln.pnext,tuples=tuples)
                        if not c2: continue
                    if contigs and not c2 in contigs: continue

                    if skipI and c1 == c2: continue

                    if ocontig2 in self.masked_regions:
                        overlapping_masked_segs = [ (aa,bb) for aa,bb in self.masked_regions[ocontig2] if aa<aln.pnext and aln.pnext<bb ]
                        if overlapping_masked_segs: 
                            if debug: print("#m2",ocontig,ocontig2, aln.pos,aln.pnext,overlapping_masked_segs ,aln.query_name,sep="\t")
                            continue  #the pair maps to a masked region

                    if not (c1,c2) in links: 
                        links[c1,c2]=[]
                    links[c1,c2].append((x1,x2))
                    if debug: print("#p",c1,c2,aln.query_name,aln.pos,aln.pnext,sep="\t")
                    #                print(self.broken_coord(ocontig,aln.pos),self.broken_coord(bam.getrname(aln.rnext),aln.pnext))
#        for c1,c2 in links:
#            print(c1,c2,links[c1,c2])
        return(links)

    def validate(self):
        contig_mins={}
        contig_maxs={}
        scaffold_mins={}
        scaffold_maxs={}
        contig_lengths={}
        counts={}
        #        scaffold_lengths={]
        lastx=0
        last_scaffold=False
        last_contig=(False,)
        for scaffold,contig,base,end,x,l in self.layout_lines:
            if last_scaffold and scaffold == last_scaffold and not x>lastx: 
                print(last_scaffold,scaffold,x,lastx)
                raise Exception 
            if x<0: raise Exception
            counts[contig,base]=counts.get((contig,base),0)+1
            contig_mins[contig,base]=min( contig_mins.get((contig,base),1e99) , x )
            contig_maxs[contig,base]=max( contig_maxs.get((contig,base),-1) , x )
            if (contig,base) in contig_lengths:
                if not contig_lengths[contig,base]==l: raise Exception 
            else:
                contig_lengths[contig,base]=l
            if (contig,base) == last_contig:
#                print(scaffold,contig,base,end,x,l,last_contig,lastx,x-lastx,"XXXX",sep="\t")
                if not l == x-lastx :
                    print(scaffold,contig,base,end,x,l,last_contig,lastx,sep="\t")
                    raise Exception 
            scaffold_mins[scaffold]=min( scaffold_mins.get(scaffold,1e99) , x )
            scaffold_maxs[scaffold]=max( scaffold_maxs.get(scaffold,-1)   , x )
            last_contig=(contig,base)
            lastx,last_scaffold=x,scaffold
        for contig,base in counts.keys():
            if not counts[contig,base]==2: raise Exception 
            if not contig_maxs[contig,base]-contig_mins[contig,base]==contig_lengths[contig,base]: 
                print("wtf?",contig,base,contig_mins[contig,base],contig_maxs[contig,base],contig_lengths[contig,base],contig_maxs[contig,base]-contig_mins[contig,base])
                raise Exception 

        self.setup_mapper()
        lastk=(1,0)
        lastl=0
        lastx=0
        for k in sorted(self.contigx.keys()):
            if k[0]==lastk[0]:
                #print("XX",k,contig_lengths.get(k),self.contigx[k],lastl,k[1]-lastk[1],lastl==k[1]-lastk[1],sep="\t")
                if not lastl==k[1]-lastk[1]:
                    print(k,contig_lengths.get(k),self.contigx[k],lastl,k[1]-lastk[1],lastl==k[1]-lastk[1],sep="\t")
                    raise Exception 
            lastk=k
            lastl = contig_lengths.get(k)
            lastx=self.contigx[k]
#        self.contigx={}
#        self.contigst={}


        print("#validation success")
#            contig_maxs[contig,base]=max( contig_maxs.get((contig,base),-1) , x )


    def add_breakpoints_ranges(self,breaks,debug=False,scores={},unlink_threshold=100.0, contig_break_threshold=100.0,score_types={},clipped_likelihood_for_gaps_only=True):
        """This function filters breaks expressed as spans according to several criteria before passing selected break points along to add_breakpoints()"""

        if not self.layout_lines or len(self.layout_lines)==0:
            self.make_trivial_layout_lines(debug=debug)
            self.validate()

        break_points = []
        category={}
        br={}
        scaffold_breaks={}
        for s,a,b,c in breaks:
            if not s in scaffold_breaks: scaffold_breaks[s]=[]

        for s,a,b,c in breaks:
            scaffold_breaks[s].append((a,b,c))

        for s in scaffold_breaks.keys():
             scaffold_breaks[s].sort()
        i=0
        offset={}
        while i< len(self.layout_lines)-1:
            #print(self.layout_lines[i])
            scaffold1,contig1,base1,end1,x1,l21 = self.layout_lines[  i  ]
            scaffold2,contig2,base2,end2,x2,l22 = self.layout_lines[ i+1 ]
            i+=1
            if not scaffold1 in scaffold_breaks: continue
            if not scaffold1 == scaffold2 : continue
            for a,b,c in scaffold_breaks[scaffold1]:
                if pairs_overlap((a,b),(x1,x2)) and (not (contig1,base1)==(contig2,base2)):  # you don't have to break a contig to make this break
                    if x1<=c and c<=x2:
                        this_offset_from_min_score = 0
                    else:
                        this_offset_from_min_score = min(abs(x1-c),abs(x2-c))

                    if (not (scaffold1,a,b,c) in category) or \
                       (category[scaffold1,a,b,c]=="break") or \
                       (category[scaffold1,a,b,c]=="gap" and offset[scaffold1,a,b,c]>this_offset_from_min_score):
                        category[scaffold1,a,b,c]="gap"
                        br[scaffold1,a,b,c]=int((x1+x2)/2)
                        offset[scaffold1,a,b,c] = this_offset_from_min_score
                elif x1<a and b<x2:
                    if not (scaffold1,a,b,c) in category:
                        if scores.get((scaffold1,a,b,c),-1e9) < contig_break_threshold:
                            category[scaffold1,a,b,c]="break"
                            br[scaffold1,a,b,c]=c
        categories={}
        for s,a,b,c in br.keys():
            if (category[s,a,b,c]=="break" and scores.get((s,a,b,c),100)< contig_break_threshold) or \
               (category[s,a,b,c]=="gap"   and scores.get((s,a,b,c),100)< unlink_threshold) :

                if clipped_likelihood_for_gaps_only and (score_types.get((s,a,b,c),"fine")=="clipped") and category[s,a,b,c]=="break": 
                    if debug: print("Skip break because we're only using clipped likelihood for un-scaffolding, not contig breaking",s,a,b,c,scores.get((s,a,b,c),-100),sep="\t")
                    if debug: print("X",s,a,b,c,br[s,a,b,c],category[s,a,b,c],score_types.get((s,a,b,c),"fine"),scores.get((s,a,b,c),-100))
                    continue
                
                break_points.append((s,br[s,a,b,c]))
                if True or debug: print("Q",s,a,b,c,br[s,a,b,c],category[s,a,b,c],score_types.get((s,a,b,c),"fine"),scores.get((s,a,b,c),-100))
                categories[s,br[s,a,b,c]]=category[s,a,b,c]
        self.add_breakpoints(break_points,debug=debug)

    def add_breakpoints(self,breaks,debug=False,categories={}):
        """breaks is a list of (scaffold,x) pairs where scaffold is the id of a scaffold in the assembly, 
and x is a basepair offset from the scaffold 5' end of the scaffold where a break should be added.  If x
falls in a gap between contigs, the linkage is broken.  If x falls in a contig, the contig is broken."""
        import networkx as nx

        g = nx.DiGraph()

        scaffold_breaks={}

        for s,x in breaks:
            if not s in scaffold_breaks: scaffold_breaks[s]=[]

        for s,x in breaks:
            scaffold_breaks[s].append(x)

        for s in scaffold_breaks.keys():
             scaffold_breaks[s].sort()

        if debug: print(scaffold_breaks)

        new_layout_lines=[]

        scaffold_broken=False

        break_indices={}

        i=0
        break_counter=0
#        debug=False

        edge_breaks={}
        elength={}
        scaffold={}
        xcoord={}
#
#  Convert the layout lines into a directed graph
#
        contig_length={}
        while i< len(self.layout_lines)-1:
            #print(self.layout_lines[i])
            scaffold1,contig1,base1,end1,x1,l21 = self.layout_lines[  i  ]
            scaffold2,contig2,base2,end2,x2,l22 = self.layout_lines[ i+1 ]
            scaffold[contig1,base1]=scaffold1
            scaffold[contig2,base2]=scaffold2
            i+=1
            if scaffold1==scaffold2: 
                if scaffold1 in scaffold_breaks and (not (contig1,base1)==(contig2,base2)) and len([z for   z in scaffold_breaks[scaffold1] if x1<z and z<x2])>0:  
                    pass # this is an OO join to remove
                else:
                    g.add_edge( (contig1,base1,end1),(contig2,base2,end2) )
                    elength[ (contig1,base1,end1),(contig2,base2,end2) ]=x2-x1
                    if (contig1,base1)==(contig2,base2): contig_length[contig1,base2]=x2-x1

                    if not ( scaffold1==scaffold2 and scaffold1 in scaffold_breaks):
                        continue

                    else:
                        break_positions = [z for   z in scaffold_breaks[scaffold1] if x1<z and z<x2]

                        if len(break_positions)==0: 
                            continue

                        edge_breaks[(contig1,base1,end1),(contig2,base2,end2)]=[]
                        for x in break_positions: # these's a break in here!

                            edge_breaks[(contig1,base1,end1),(contig2,base2,end2)].append( x-x1 )

        for e1,e2 in edge_breaks.keys():
#            print(e1,e2)
            next_nodes= g.neighbors(e2) 
            prev_nodes= g.predecessors(e1) 
#            print(prev_nodes)
            if len(next_nodes)>1: 
                raise Exception 
            if len(prev_nodes)>1: 
                raise Exception 
            next_node=False
            prev_node=False
            if next_nodes: next_node = next_nodes[0]
            if prev_nodes: prev_node = prev_nodes[0]
            g.remove_edge(e1,e2)

            #print(e1,e2,edge_breaks[e1,e2],scaffold[e1[0],e1[1]],next_node,prev_node)

            if not (e1[0],e1[1])==(e2[0],e2[1]):
                #print("this is just an oo link to udo")
                #this is just an edge to remove
                continue

            ocontig = e1[0]
            
            if e1[2] in ["5",5]: #forward

                
                
                if next_node: g.remove_edge(e2,next_node)
                g.remove_node(e2)

                #                          a        b
                #         prev    e1     x1          x2       e2   next
                #  ---------->     --------------------------->    ----------------->

                #just implement the first and last break.  (this is not necessarily good if the contig is really, really long. XXX TODO )
                edge_breaks[e1,e2].sort()
                x1=edge_breaks[e1,e2][0]
                x2=edge_breaks[e1,e2][-1]

                for i in range(len(edge_breaks[e1,e2])-1):
                    xa = edge_breaks[e1,e2][i  ]
                    xb = edge_breaks[e1,e2][i+1]
                    new_node_a = ( e1[0] , e1[1]+xa, '5' )
                    new_node_b = ( e1[0] , e1[1]+xa, '3' )
                    g.add_edge(new_node_a,new_node_b)
                    elength[ new_node_a,new_node_b ] = xb-xa
                    scaffold[e1[0] , e1[1]+xa] = scaffold[e1[0],e1[1]]

                new_node_1 = (e1[0],e1[1]   ,'3')
                new_node_2 = (e1[0],e1[1]+x2,'5')
                new_node_3 = (e1[0],e1[1]+x2,'3')
                scaffold[e1[0],e1[1]+x2] = scaffold[e1[0],e1[1]]

                if debug: print("b",e1,e2,prev_node,next_node,x1,x2,edge_breaks[e1,e2],new_node_1,new_node_2,   sep="\t")

                g.add_edge(e1,new_node_1)
                g.add_edge(new_node_2,new_node_3)
                if next_node: g.add_edge(new_node_3,next_node)

                LL = elength[e1,e2]
                if next_node: LL2 = elength[e2,next_node]
                elength[ e1,new_node_1 ] = x1
                elength[ new_node_2,new_node_3 ] = LL - x2
                if next_node: elength[ new_node_3,next_node  ] = LL2


                if debug: print("b",e1,e2,prev_node,next_node,x1,x2,edge_breaks[e1,e2],new_node_1,new_node_2, elength.get(( e1,new_node_1 )),elength.get(( new_node_2,new_node_3 )), elength.get((new_node_3,next_node)), sep="\t")
#                print("zzz",contig_length[e1[0],e1[1]],LL,elength[e1,e2],edge_breaks[e1,e2],x1,x2,elength[e1,e2],elength[ e1,new_node_1 ],elength[ new_node_2,new_node_3 ])

            else:                #reverse
                
                if prev_node: g.remove_edge(prev_node,e1)
                g.remove_node(e1)

                #just implement the first and last break.  (this is not necessarily good if the contig is really, really long. XXX TODO )
                edge_breaks[e1,e2].sort()                
                x1=elength[e1,e2]-edge_breaks[e1,e2][-1]
                x2=elength[e1,e2]-edge_breaks[e1,e2][0]


                #                          b        a
                #         prev    e1     x2          x1       e2   next
                #  ---------->    <----------------------------    ----------------->
                #                  A      B           C

                for i in range(len(edge_breaks[e1,e2])-1):
                    xb = elength[e1,e2]-edge_breaks[e1,e2][i  ]
                    xa = elength[e1,e2]-edge_breaks[e1,e2][i+1]
                    new_node_a = ( e1[0] , e2[1]+xa, '5' )
                    new_node_b = ( e1[0] , e2[1]+xa, '3' )
                    g.add_edge(new_node_a,new_node_b)
                    elength[ new_node_a,new_node_b ] = abs(xb-xa)
                    scaffold[e1[0] , e2[1]+xa] = scaffold[e1[0],e2[1]]

                new_node_C = (e1[0],e2[1]   ,'3') #
                new_node_B = (e1[0],e2[1]+x2,'5') #
                new_node_A = (e1[0],e2[1]+x2,'3') #
                scaffold[e1[0],e2[1]+x2] = scaffold[e1[0],e2[1]]

                g.add_edge(new_node_C,e2) #
                g.add_edge(new_node_A,new_node_B) #
                if prev_node: g.add_edge(prev_node,new_node_A) #

                LL =  elength[e1,e2]
                if prev_node: LL2 = elength[prev_node,e1]
                elength[ new_node_C, e2 ] = x1
                elength[ new_node_A,new_node_B ] = LL - x2
                if prev_node: elength[ prev_node, new_node_A  ] = LL2
#                print("zzW",edge_breaks[e1,e2],x1,x2,elength[e1,e2],elength[ new_node_1, e2 ],elength[ new_node_3,new_node_2 ])
            
#        for e1,e2 in elength.keys():
#            if elength[e1,e2]<0:
#                print("wtf?",e1,e2,elength[e1,e2])
#
#  Compute the new layout coordinates and scaffold IDs
#
        scaffold_incr={}
        node_x={}
        contig_length={}
        new_scaffold={}
        for n in sorted(g.nodes()):
            #print(n,scaffold[n[0],n[1]],g.degree(n),g.in_degree(n),g.out_degree(n))
            if g.in_degree(n)==0:
                scaffold_id=scaffold[n[0],n[1]]
                new_scaffold_id = scaffold_id 
                if scaffold_id in scaffold_incr: 
                    new_scaffold_id = scaffold_id + "."+str(scaffold_incr.get(scaffold_id,0))
                new_scaffold[n]=new_scaffold_id
                scaffold_incr[scaffold_id] = scaffold_incr.get(scaffold_id,0)+1
                ll=list(nx.dfs_preorder_nodes(g,n))
                xx=0
                lastl=False
                for l in ll:
                    if lastl:
                        xx+=elength[ lastl,l ]
                        #print("z",l,lastl,elength[lastl,l],node_x[lastl],xx,new_scaffold[n],sep="\t")
                        if (lastl[0],lastl[1])==(l[0],l[1]):
                            contig_length[lastl[0],lastl[1]]=elength[ lastl,l ]
                    #print(new_scaffold_id,l,xx)
                    node_x[l]=xx
                    lastl = l
#
#  Build the new layout_lines list
#
        self.scaffold_lengths={}
        new_layout = []
        for n in g.nodes():
            if g.in_degree(n)==0:

                new_scaffold_id = new_scaffold[n]
                ll=list(nx.dfs_preorder_nodes(g,n))

                for l in ll:
#                    print("w",new_scaffold_id,l[0],l[1],l[2],node_x[l])
#                    print("z",new_scaffold_id,l[0],l[1],l[2],node_x[l],contig_length[l[0],l[1]])
                    new_layout.append( (new_scaffold_id,l[0],l[1],l[2],node_x[l],contig_length[l[0],l[1]] ) )
#                    new_layout.append( (new_scaffold_id,l[0],l[1],l[2],node_x[l],elength[l[0],l[1]] ) )
                    self.scaffold_lengths[new_scaffold_id] = max( node_x[l] , self.scaffold_lengths.get(new_scaffold_id,0))

                    
        
        self.layout_lines=new_layout
#        self.layout_lines.sort( key=lambda x: (x[0],x[4]) )
        self.index_ocontigs()

    def load_assembly(self,infile):
        self.scaffold_lengths={}
        self.bams=[]
        self.layout_lines=[]
        self.model_params={}
        self.ocontig_lengths={}
        self.contig_lengths={}
        self.scaffold_ocontigs={}
        self.bam_objects=False
#        self.model_params={}
        f=open(infile,"r")
        path = ""
        #if infile.startswith("/"):
        path = os.path.dirname(os.path.normpath(infile))
        for l in f:
            #print(l)
            if l[0]=="M":
                c=l.strip().split()
                self.model_params[c[1]] = eval(" ".join(c[2:]))[0]
            if l[0]=="B":
                #self.bams.append(l[2:].strip())
                possible_path = l[2:].strip()
                if possible_path.startswith("/"):
                   self.bams.append(l[2:].strip())
                else:
                    self.bams.append(os.path.join(path,l[2:].strip()))   
            if l[0]=="S":
                #self.shotgun.append(l[2:].strip())
                possible_path = l[2:].strip()
                if possible_path.startswith("/"):
                    self.shotgun.append(l[2:].strip())
                else:
                    self.shotgun.append(os.path.join(path,l[2:].strip()))
            if l[0]=="P": # Layout
                self.layout_lines.append( parse_pline2(l,self.scaffold_lengths) )
                scaffold,contig,base,end,x,lll = self.layout_lines[-1]
                self.contig_lengths[contig+"_"+str(base)]=lll
                self.contig_lengths[contig,base]=lll
                self.scaffold_ocontigs[scaffold] = self.scaffold_ocontigs.get(scaffold,[])+[contig]
            if l[0]=="D": # Masked regions (from D for shotgun depth)
                dummy,ocontig,a,b = l.strip().split()
                if not ocontig in self.masked_regions: self.masked_regions[ocontig]=[]
                self.masked_regions[ocontig].append((int(a),int(b)))
            if l[0]=="L": # Length of input contigs, same as bam header lengths  
                dummy,ocontig,l = l.strip().split()
                #if not ocontig in self.masked_regions: self.masked_regions[ocontig]=[]
                self.ocontig_lengths[ocontig] = int(l)
        f.close()
        self.index_ocontigs()
#        print("setup mapper")
        self.setup_mapper(debug=True)

    def save_assembly(self,outfile):
        f=open(outfile,"w")
        for b in self.shotgun:
            f.write("S {}\n".format(b))
        for b in self.bams:
            f.write("B {}\n".format(b))
        #f.write("M {}\n".format(self.model_params))
        if self.model_params:
            for k in sorted(self.model_params.keys()):
                f.write("M {} {}\n".format(k,(self.model_params[k],)))
        for l in self.layout_lines:
            f.write("P "+" ".join(map(str,l))+"\n")
        for oc in sorted(self.masked_regions.keys()):
            for a,b in self.masked_regions[oc]:
                f.write("D {} {} {}\n".format(oc,a,b))
        for oc in sorted(self.ocontig_lengths.keys()):
            f.write("L {} {}\n".format(oc,self.ocontig_lengths[oc]))

#        self.masked_regions={}
        f.close()

    def scaffold_mapping():
        for i in range(len(self.layout_lines)):
            scaffold1,contig1,base1,end1,x1,l11 = self.layout_lines[i  ]
            if end1 == "3":
                yield("{}_{}".format(contig1,base1),scaffold1)
        
    def as_edges(self):
        for i in range(len(self.layout_lines)-1):
            scaffold1,contig1,base1,end1,x1,l21 = self.layout_lines[i  ]
            scaffold2,contig2,base2,end2,x2,l22 = self.layout_lines[i+1]
            if scaffold1==scaffold2 :
                if contig1==contig2 and base1==base2: 
                    yield(["#edge:","{}_{}.{}".format(contig1,base1,end1) ,"{}_{}.{}".format(contig2,base2,end2) ,{'length':l21,'contig':True}])
                else:
                    yield(["#edge:","{}_{}.{}".format(contig1,base1,end1) ,"{}_{}.{}".format(contig2,base2,end2) ,{'length':x2-x1,'contig':False}])
#                    yield(["#edge:",contig1+"."+end1,contig2+"."+end2,{'length':x2-x1,'contig':False}])

    def contig_coords(self):
        seen = {}
        for l in self.layout_lines:
            scaffold1,contig1,base1,end1,x1,l21 = l
            if not (contig1,base1) in seen:
                yield( {"contig":contig1,"base":base1,"scaffold":scaffold1,"span":(x1,x1+l21),"flipped":end1 in [3,"3"] })

            seen[contig1,base1]=True

    def dump_layout(self,f):
#        f=open(outfile,"w")
        for s in self.scaffold_lengths.keys():
            f.write("{} slen\n".format(self.scaffold_lengths[s]))
        for l in self.layout_lines:
            scaffold1,contig1,base1,end1,x1,l21 = l
            f.write("p: {} {}_{}.{} {} - -1 {} {} {}\n".format(scaffold1,contig1,base1,end1,x1,self.scaffold_lengths[scaffold1],l21,False))
#        f.close()
        #p: 1 Scaffold105707_1.3 3019 - -1 85132546 3048 False

    def load_playout(self,layoutfile):
        #print(layoutfile)
        self.scaffold_lengths={}
#        self.contig_lengths={}
        if not os.path.exists(layoutfile):
            raise Exception 
        self.layout_lines = [ parse_pline(c,self.scaffold_lengths) for c in open(layoutfile,"rt") if c[:2]=="p:" ]
        #for l in self.layout_lines:
            #print(l)
        self.index_ocontigs()

    def check_bams(self):
        for bamfile in self.bams:
            if not os.path.exists(bamfile):
                raise Exception 

    def set_modelparams(self,datamodel):
        if os.path.exists(datamodel):
            self.model_params = eval(open(datamodel).read())
        else:
            self.model_params=datamodel

if __name__=="__main__":

    import sys
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d','--debug',default=False,action="store_true",help="Turn on debugging ouput")
    parser.add_argument('-E','--edges',default=False,action="store_true",help="Print the layout as edges.")
#    parser.add_argument('-E','--edges',default=False,action="store_true",help="Print the layout as edges.")
    parser.add_argument('-D','--dump_layout',default=False,action="store_true",help="Print the layout as edges.")
    #parser.add_argument('-p','--progress',default=False,action="store_true",help="Print progress info")
    parser.add_argument('-M','--datamodel',default=False, help="Name of a file containing the model parameter dictionary")
    parser.add_argument('-m','--mask'   ,default=[], action="append",help="Name of a file containing regions of the input assembly to mask out")
    parser.add_argument('-b','--bamfile',default=[], action="append",help="A chicago bam file")
    parser.add_argument('-S','--shotgunbamfile',default=[], action="append",help="A shotgun bam file")
    #parser.add_argument('-p','--progress',default=False,action="store_true",help="Print progress info")
    parser.add_argument('-L','--layout',default=False,help="File containing a scaffolding layout in 'p:' lines.")
    parser.add_argument('-o','--outfile',default=False,help="Filename for serialization of the assembly file.")
    parser.add_argument('-i','--infile',default=False,help="Filename for serialised assembly input file.")

    args = parser.parse_args()
    
#    print(args)
    if args.infile:
        asf = HiriseAssembly()
        asf.load_assembly(args.infile)
#        asf.validate()
    else:
        asf = HiriseAssembly( {"bams":args.bamfile, "datamodel": args.datamodel, "layout": args.layout , "debug":args.debug, "shotgun": args.shotgunbamfile} )
    if args.debug: print("# done with intital load")
    asf.validate()
#    asf.index_ocontigs()

    if args.datamodel:
        asf.set_modelparams(args.datamodel)

    for segments_file in args.mask:
        asf.add_mask_regions(filename=segments_file)
        asf.merge_masked_regions()

    if args.outfile:
        asf.save_assembly( args.outfile )
        
    if args.edges:
        for e in asf.as_edges():
            print("\t".join(map(str,e)))

#    asf.add_breakpoints( [( '1', 1692 ),('1',10000)] )
#    asf.add_breakpoints( [( '1', 1692 ),('1',2000),('1',10000)] )
#    asf.add_breakpoints( [( '1', 1692 ),('1',2000),('1',2100),('1',10000)] )

#    asf.validate()

#    asf.load_ocontig_info()
#    print("#loaded ocontig info")
#    asf.load_ocontig_info()
#    print("#loaded ocontig info")
#    print(asf.get_links(["Scaffold116325"],links={}))
#    print(asf.get_links(["Scaffold116325"],links={}))
#    asf.get_scaffold_links("1000",skipI=True)

#1691
    if args.dump_layout:
        asf.dump_layout(sys.stdout)

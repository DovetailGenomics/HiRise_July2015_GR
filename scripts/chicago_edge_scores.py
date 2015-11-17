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
import math
import sys
"""Chicago likelihood model calculations as described in ... """

debug=False

def fragment_size(l1,l2,o1,o2,coords,g=0):
    """o1 and o2 indicate the orientations of the two contigs being used  0 means 5' to 3'; 1 means flipped"""
    x,y = coords[0],coords[1]
    if (o1,o2) == (0,0):     #         -------1----->    ----2------>
        return (y+l1-x+g)
    if (o1,o2) == (0,1):     #         -------1----->   <---2-------
        return (l2-y+l1-x+g)
    if (o1,o2) == (1,0):     #        <-------1-----     ----2------>
        return (x+y+g)
    if (o1,o2) == (1,1):     #        <-------1----->   <----2------
        return (l2-y+x+g)

class ChicagoModel(object):
    def __init__(self,params):
        self.alpha = params.get('alpha',[1.0])
        self.beta  = params.get('beta',[1.0])
        self.pn    = params.get('pn',0.5)
        self.G     = params.get('G',3000000000)
        self.N     = params.get('N',100000000)
        self.nexps = len(self.alpha)

    def set_params(self,params):
        self.alpha = params.get('alpha')
        self.beta  = params.get('beta' )
        self.pn    = params.get('pn'   )
        self.G     = params.get('G'    )
        self.N     = params.get('N'    )
        self.nexps = len(self.alpha)
        self.lnFinf = math.log(self.pn)-math.log(self.G)
        self.logN = math.log(self.N)
        self.logG = math.log(self.G)

    def __repr__(self):
        return(str({ 'pn': self.pn, 'G':self.G, 'N':self.N, 'a':self.alpha, 'b':self.beta, 'n':self.nexps }))

### These are the things that depend on the model:
#
    def f0(self,x):
        """P_{raw}(x)= \sum \alpha_i \beta_i e^{-x \beta_i}  """
        return sum( [ self.alpha[i]*self.beta[i]*math.exp(-x*self.beta[i]) for i in range(self.nexps) ] )

    def f0_prime(self,x):
        alpha=self.alpha
        beta=self.beta
        """P_{raw}(x)= \sum \alpha_i \beta_i e^{-x \beta_i}  """
        return ( sum([ -alpha[i]*beta[i]*beta[i]*math.exp(-beta[i]*x) for i in range(self.nexps) ]) )

    def f0_double_prime(self,x):
        alpha=self.alpha
        beta=self.beta
        """P_{raw}(x)= \sum \alpha_i \beta_i e^{-x \beta_i}  """
        return  ( sum([ alpha[i]*(beta[i]**3.0)*math.exp(-beta[i]*x) for i in range(self.nexps) ]) )

    def F0(self,d):
        # \int f(x) dx
        alpha=self.alpha
        beta =self.beta
        if d<0:
            print("WTF? d<0:",d,file=sys.stderr)
            raise Exception 
        return self.pn*d/self.G + (1.0-self.pn) * sum( [ - alpha[i]*(math.exp(-beta[i]*d))  for i in range(self.nexps)] )

    def H0(self,d):
        #   \int b x f(x) dx = \int b x e^{- b x} dx = -e^{- b x} \frac{bx+1}{b}
        alpha=self.alpha
        beta =self.beta
        return  0.5*d*d*self.pn/self.G + (1.0-self.pn)* sum([ -(old_div(alpha[i],beta[i]))*(math.exp(-beta[i]*d)*(beta[i]*d + 1.0)) for i in range(self.nexps) ])

    def H0_good(self,d):
        #   \int b x f(x) dx = \int b x e^{- b x} dx = -e^{- b x} \frac{bx+1}{b}
        alpha=self.alpha
        beta =self.beta
        return  (1.0-self.pn)* sum([ -(old_div(alpha[i],beta[i]))*(math.exp(-beta[i]*d)*(beta[i]*d + 1.0)) for i in range(self.nexps) ])

#
###

    def F(self,d):
        return self.F0(d)-self.F0(0)

    def H(self,d):
        return self.H0(d)-self.H0(0)

    def f(self,x,cache={}):   #cache this
        if x in cache: return cache[x]
        r= old_div(self.pn,self.G) + (1.0-self.pn) * self.f0(x)
        cache[x]=r
        return r

    def f_prime(self,x):
        return (1-self.pn)*self.f0_prime(x)

    def f_double_prime(self,x):
        return (1-self.pn)*self.f0_double_prime(x)

### These are for backwards compatability with the old implementation
#
    def p_insert_raw(self,x):
        return self.f0(x)

    def p_insert_effective(self,x):
        """ p_n/G + (1-p_n) P_{raw}(x) """
        return  self.f(x)

    def ddg_p(self,x):
        return self.f0_prime(x)

    def d2dg2_p(self,x):                                                                                                                                                     
        return self.f0_double_prime(x)

    def ll(self,l1,l2,o1,o2,links,g,p0=-1):
        return self.lnL(l1,l2,o1,o2,g,links)

    def p0(self,l1,l2,g):
        return 1.0-(old_div(self.n_bar(l1,l2,g),self.N))
#
###

    def omega(self,l1a,l2a,d):
        l1=min(l1a,l2a)
        l2=max(l1a,l2a)
        if     d<0:  return 1
        elif   d<l1: return d+1
        elif   d<l2: return l1+1
        elif   d<=l1+l2: return l1+l2-d+1
        else: return 0

    def osigma(self,l1a,l2a,d):
        l1=min(l1a,l2a)
        l2=max(l1a,l2a)
        if     d<0:  return 0
        elif   d<l1: return 1.0
        elif   d<l2: return 0
        elif   d<l1+l2: return -1.0
        else: return 0

    def lnF(self,x,cache={}):   # Cache this
        if x in cache: return cache[x]
        r=math.log( self.f(x) )
        cache[x]=r
        return r  #math.log( self.f(x) )

    def T(self,d,cache={}):  # Cache this
        if d in cache: return cache[d]
        #\sum_{i=0}^{d} (d-i) f(i) \approx 
        x= d*self.F(d) - self.H(d)
        cache[d]=x
        return x

    def p(self,l1,l2,g):
        if l1<0: raise Exception 
        if l2<0: raise Exception 
        if (l1+l2+g)<0: sys.stderr.write("wtf: l1+l2+g < 0 ; {}+{}+{} < 0\n".format(l1,l2,g))
        if (l1+g)<0: sys.stderr.write("wtf: l1+g < 0 ; {}+{} < 0\n".format(l1,g))
        if (l2+g)<0: sys.stderr.write("wtf: l2+g < 0 ; {}+{} < 0\n".format(l2,g))
        if (g)<0: sys.stderr.write("wtf: g < 0 ; {} < 0\n".format(g))
        p = old_div((self.T(l1+l2+g)+self.T(g)-self.T(l1+g)-self.T(l2+g)),self.G)
        return p
        
    def p_prime(self,l1,l2,g):
        if l1<0: raise Exception 
        if l2<0: raise Exception 
        p = old_div((self.F(l1+l2+g)+self.F(g)-self.F(l1+g)-self.F(l2+g)),self.G)
#        p = ((l1+l2+g)*self.f(l1+l2+g)+g*self.f(g)-(l1+g)*self.f(l1+g)-(l2+g)*self.f(l2+g))/self.G
        return p

    def p_double_prime(self,l1,l2,g):
        if l1<0: raise Exception 
        if l2<0: raise Exception 
        ss =  self.f(l1+l2+g)+self.f(g)-self.f(l1+g)-self.f(l2+g)
#        ss +=(l1+l2+g)*self.f_prime(l1+l2+g)+g*self.f_prime(g)-(l1+g)*self.f_prime(l1+g)-(l2+g)*self.f_prime(l2+g)
        return old_div(ss,self.G)

#    def Q(self,l1,l2,g,x):
#        return (self.omega(l1,l2,x-g)*self.f_prime(x) + self.osigma(l1,l2,x-g)*self.f(x) )/( self.omega(l1,l2,x-g)*self.f(x) )

    def R(self,l1,l2,g,x):
        # P''/P - (P'/P)^2
        return ( old_div(self.f_double_prime(x),self.f(x))  - (self.Q(x))**2.0 )

    def n_bar0(self,l1l2):
        #return self.N*self.pn*l1l2*2/(self.G*self.G)
#        return self.N*self.pn*l1l2*2.0/(self.G*self.G)
        return self.N*self.pn*l1l2/(self.G*self.G)

    def n_bar(self,l1,l2,g):
        p = self.p(l1,l2,g)
        return self.N*p

    def cutScore(self,l1,x,n,sumLogP,rangeCutoff=False,minSep=0):
        la=x-minSep/2
        lb=l1-x-minSep/2
        if rangeCutoff and la>rangeCutoff: la=rangeCutoff
        if rangeCutoff and lb>rangeCutoff: lb=rangeCutoff
        if la <=0: return 0.0
        if lb <=0: return 0.0
        try:
            n_bar = self.n_bar(la,lb,minSep)
        except Exception as e:
            print("Caught an exception when trying to compute n_bar, where la={}, lb={} and gap={}, based on l1={}, x={}, minSep={}".format(la,lb,minSep,l1,x,minSep))
            raise e
        n_bar0= self.n_bar0(la*lb)
        sumLogP0 = n*(math.log(self.pn) - self.logG)
        return -n_bar + sumLogP - ( -n_bar0 + sumLogP0)

    def tileScore(self,binwidth,tile,n,sumLogP,rangeCutoff=False,minSep=0):
        la=binwidth
        lb=binwidth
#        if rangeCutoff and la>rangeCutoff: la=rangeCutoff
#        if rangeCutoff and lb>rangeCutoff: lb=rangeCutoff
        if tile[0]==tile[1]:
            n_bar = self.N*(binwidth*(self.F(binwidth)-self.F(minSep))-(self.H(binwidth)-self.H(minSep)))/self.G
            n_bar0= self.N*self.pn*(((binwidth-minSep)/self.G)**2.0)
#            sumLogP0 = n*(math.log(self.pn) - self.logG)
            sumLogP0 = n*self.lnFinf
            return -n_bar + sumLogP - ( -n_bar0 + sumLogP0)        
            
        else:
            n_bar = self.n_bar(la,lb,(tile[1]-tile[0]-1)*binwidth)
            n_bar0= self.N*self.pn*la*lb*2/(self.G*self.G)
#            sumLogP0 = n*(math.log(self.pn) - self.logG)
#        self.lnFinf = math.log(self.pn)-math.log(self.G)
            sumLogP0 = n*self.lnFinf
            return -n_bar + sumLogP - ( -n_bar0 + sumLogP0)        

    def lnL(self,l1,l2,o1,o2,g,links):
        n=len(links)
        n_bar = self.n_bar(l1,l2,g)
        try:
            r= n*self.logN - n_bar - n*self.logG + sum( [ math.log(self.omega(l1,l2,fragment_size(l1,l2,o1,o2,links[i],0))) + self.lnF( fragment_size(l1,l2,o1,o2,links[i],g) ) for i in range(n) ] )
        except Exception as e:
            print(e)
            print("l1,l2",l1,l2)
            print("o1,o2",o1,o2)
            print("links:",links)
            print("fragment sizes g=0:", [fragment_size(l1,l2,o1,o2,links[i],0) for i in range(n)])
            print("fragment sizes g=g:", [fragment_size(l1,l2,o1,o2,links[i],g) for i in range(n)])
            print("omega=",[self.omega(l1,l2,fragment_size(l1,l2,o1,o2,links[i],0)) for i in range(n)])
            print("lnF=",[ self.lnF( fragment_size(l1,l2,o1,o2,links[i],g)) for i in range(n)])
            raise e
        return r

    def lnL0(self,l1,l2,o1,o2,links):
        n=len(links)
        n_bar = self.n_bar0(l1*l2)
#        n_bar = self.N*self.pn*l2*l1/(self.G**2)
#        return n*math.log(n_bar) - n_bar - math.log(math.factorial(n))
        ll1= n*self.logN - n_bar - n*self.logG + sum( [ math.log(self.omega(l1,l2,fragment_size(l1,l2,o1,o2,links[i],0))) + self.lnFinf for i in range(n) ] )  
        return ll1
#        ll2=ll1
#        o2+=1
#        o2=o2%2
#        ll2= n*self.logN - n_bar - n*self.logG + sum( [ math.log(self.omega(l1,l2,fragment_size(l1,l2,o1,o2,links[i],0))) + self.lnFinf for i in range(n) ] )  
#        return max(ll1,ll2)
 
    def score(self,l1,l2,o1,o2,links,g,p0=0):
        n=len(links)
        #return self.lnL(l1,l2,o1,o2,g,links) - self.lnL0(l1,l2,o1,o2,links)
        # alternatively:
        return( - self.n_bar(l1,l2,g) + self.n_bar0(l1*l2) - n * self.lnFinf + sum( [ ( self.lnF( fragment_size(l1,l2,o1,o2,links[i],g) ) ) for i in range(n) ] )  )

    def Q(self,x):
#        return (self.omega(l1,l2,x-g)*self.f_prime(x) + self.osigma(l1,l2,x-g)*self.f(x) )/( self.omega(l1,l2,x-g)*self.f(x) )
        return (old_div(self.f_prime(x),self.f(x)))

    def ddg_lnL(  self,l1,l2,o1,o2,links,g):
        n=len(links)

        try:
            r =  -self.N * self.p_prime(l1,l2,g) + sum([ self.Q(fragment_size(l1,l2,o1,o2,links[i],g) ) for i in range(n)])
        except Exception as e:
            print(e)
            print("l1,l2,g",l1,l2,g)
            print("n_bar:",self.n_bar(l1,l2,g))
            print("p:",self.p(l1,l2,g))
            print("p_prime:",self.p_prime(l1,l2,g)) 
            raise e
        return r

#        return ( -self.n_bar(l1,l2,g) * self.p_prime(l1,l2,g) / self.p(l1,l2,g) + sum([ self.Q(l1,l2,g,fragment_size(l1,l2,o1,o2,links[i],g) ) for i in range(n)]))

    def d2dg2_lnL(self,l1,l2,o1,o2,links,g):
        n=len(links)
        return ( -self.N * self.p_double_prime(l1,l2,g) + sum([ self.R(l1,l2,g,fragment_size(l1,l2,o1,o2,links[i],g) ) for i in range(n) ]) )

    def d2dg2_llr(self, l1,l2,o1,o2,links,g,seen={}):
        if not "warned" in seen: 
            sys.stderr.write("using deprecated interface to likelihood model code.\n")
            seen["warned"]=1
        return self.d2dg2_lnL(l1,l2,o1,o2,links,g)

    def ddg_llr(self, l1,l2,o1,o2,links,g,seen={}):
        if not "warned" in seen: 
            sys.stderr.write("using deprecated interface to likelihood model code.\n")
            seen["warned"]=1
        return self.ddg_lnL(l1,l2,o1,o2,links,g)

    def ml_gap(self,l1,l2,o1,o2,links,g0):
        gap=g0
        G=self.G
        N=self.N
        pn=self.pn

        last_gap=g0
        for i in range(100):
            #p0 = self.p0(      l1,l2,gap)
            x1=  self.ddg_lnL(  l1,l2,o1,o2,links,gap) #1st derivative of the likelihood wrt. gap size
            x2=  self.d2dg2_lnL(l1,l2,o1,o2,links,gap) #2nd derivative of the likelihood wrt. gap size
            score=self.score(l1,l2,o1,o2,links,gap)
#            if debug: print "\t".join(map(str,[ "#it",i,o1,o2,gap,score,gap-x1/2,x1,x2,x1/x2]))
            if x2==0:                     
                break                     
                print("hit x2==0 after",i)  
            gap = int( gap - old_div(x1,x2) )
            if gap<0.0:
                gap=10.0

            if gap>200000.0:
                gap=200000

            if abs(old_div(x1,x2))<0.1: break
            if abs(gap-last_gap)<1.0: break
            last_gap=gap

        score=self.score(l1,l2,o1,o2,links,gap) #l1,l2,o1,o2,G,pn,links,N,gap,p0)
        
        if gap<0.0:
            print("wtf? negative gap")
            raise Exception 

        return gap,score
    
model=ChicagoModel({})

def insert_size_dist(x):
    return model.p_insert_raw(x)

def set_exp_insert_size_dist_fit_params(fit_params):
    model.set_params(fit_params)
    #print("#",model)

def p_not_a_hit(l1,l2,GenomeSize,gaplen,pn):
    return model.p0(l1,l2,gaplen)

def llr_v0(l1,l2,o1,o2,GenomeSize,pn,links,N,gaplen,p0 ):
    return model.score(l1,l2,o1,o2,links,gaplen,p0)

def ll(l1,l2,o1,o2,GenomeSize,pn,links,N,gaplen,p0=-1):
    return model.ll(l1,l2,o1,o2,links,gaplen,p0)

def ml_gap(l1,l2,o1,o2,G,pn,links,N,g0):
    #print "#",links
    return model.ml_gap(l1,l2,o1,o2,links,g0)
#   def ml_gap(self,    l1,l2,o1,o2,links,g0):



if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-M','--set_insert_size_dist_fit_params',default=False )
    parser.add_argument('-d','--debug',default=False )

    args = parser.parse_args()
    debug=args.debug

    fmodel=open( args.set_insert_size_dist_fit_params )
    contents = fmodel.read()
    try:
        fit_params=eval(contents)
    except:
        "couldn't deal with option", args.param
    fmodel.close
    set_exp_insert_size_dist_fit_params(fit_params)

    for x in range(1,200000,100):
        print(x,model.f(x),model.f_prime(x),model.S(x),model.F(x),model.H(x))
    


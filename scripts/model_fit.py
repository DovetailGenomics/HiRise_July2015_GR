#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from past.utils import old_div
import sys
import math
import random
import json

def log(s):
    sys.stderr.write(  "%s\n" %(s) )

def fx(x,model):
    a=model['a']
    b=model['b']
    nu=model['nu']
    l=model['l']
    n=len(a)
    if not len(b)==n: 
        raise Exception 
    p= ( nu + l * sum( [ a[i]*b[i]*math.exp(-b[i]*(x)) for i in range(n)] ) )
#    if p < 1.0e-9: p = 1.0e-9
    return ( p ) 


def model_string(Ne,model,pn,G,n,params):
    fparams=dict(params)
    fparams['N']=Ne
    fparams['beta']=model['b']
    fparams['alpha']=model['a']
    fparams['pn']=pn
    fparams['fn']="{}*({}+(1.0-{})*(".format(Ne,old_div(pn,G),pn)+ " + ".join(["{}*exp(-(x)/{})".format( model['a'][i]*model['b'][i],old_div(1.0,model['b'][i]) ) for i in range(n) ])+"))"
    return(str(fparams))

def model_json(Ne,model,pn,G,n,params):
    fparams=dict(params)
    fparams['N']=Ne
    fparams['beta']=model['b']
    fparams['alpha']=model['a']
    fparams['pn']=pn
    fparams['fn']="{}*({}+(1.0-{})*(".format(Ne,old_div(pn,G),pn)+ " + ".join(["{}*exp(-(x)/{})".format( model['a'][i]*model['b'][i],old_div(1.0,model['b'][i]) ) for i in range(n) ])+"))"
    return(fparams)

def mutate(model,stepsize):
    model2=dict(model)

    if random.random()<0.8:
        n=len(model['a'])
        n2=[(random.random()-0.5)*2.0*stepsize for i in range(n) ]
        n3 = [ max(0.0,n2[i] + model['a'][i]) for i in range(n) ]
        a=sum(n3)
        model2['a'] = [old_div(n3[i],a) for i in range(n)] 
    else:
        model2['l']+=(random.random()-0.5)*5.0e6
        model2['l']=max(0.0,model2['l'])
    return model2

    
def chisq(model,data,xmin,xmax=1e6):
    chi=0
    for x,y in list(data.items()):
        if x<xmin: continue
        if x>xmax: continue

#        print "#",x,y,model
        z=fx(x,model)
        if y<=0.0: 
            print("#wtf y",y)
            raise Exception 
        if z<=0.0: 
            print("#wtf z",z,x,y)
            print(model)
            raise Exception 
        chi+=(math.log(y)-math.log(z))**2.0
    return chi

if __name__=="__main__":

    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-o','--outfile')
    parser.add_argument('-N','--nnoise',required=False,type=float,help="estimated number of noise reads; overrides contents of -P 'Nn'")
    parser.add_argument('-G','--genomesize',type=float,default=3.0e9,help="Genome size; overrides contents of -P 'G'")
    parser.add_argument('-s','--steps',required=False,type=int,default=3000)
    parser.add_argument('-C','--nexponential_components',required=False,type=int,default=10)
    parser.add_argument('-m','--min',required=False,type=float,default=1000.0)
    parser.add_argument('-x','--max',required=False,type=float,default=200000.0)
    parser.add_argument('-P','--param',required=False)
    parser.add_argument('-d','--debug',default=False,action="store_true")
    parser.add_argument('-p','--progress',default=False,action="store_true")
    parser.add_argument('--seed',required=False,type=int,default=1, help="Seed for random number generation, use -1 for no seed")
    parser.add_argument('--json', action="store_true", default="False")
    args = parser.parse_args()
    if args.seed != -1 :
      random.seed(args.seed)

    fit_params={}
    try:
        fit_params = eval(args.param)
    except Exception as e:
        f=open(args.param)
        contents = f.read()

        try:
            fit_params=eval(contents)
        except:
            "couldn't deal with option", args.param
            #exit(0)
        f.close


    f={}

    G=fit_params.get('G',3.0e9)
    Nn=fit_params.get('Nn',0.0)
    if args.nnoise:
        Nn=args.nnoise
    if args.genomesize:
        G=args.genomesize

    while True:
        l=sys.stdin.readline()
        if not l:break
        c=l.strip().split()
        f[int(float(c[0]))]=float(c[1])

    nu=old_div(Nn,G)

    n=args.nexponential_components
    xmin=args.min
    incr = ( args.max / xmin)**(old_div(1.0,(n-1)))
    b = tuple([ old_div(2.0,(xmin * (incr**i)))  for i in range(n) ]  )
    a = tuple([  old_div(1.0,(-math.log(b[i])))        for i in range(n)] )

    model={ 'nu':nu, 'l':Nn, 'a':a, 'b':b, 'n':n }

#    maxstep=100
    step=0
    chi = chisq(model,f,xmin)
    while step<args.steps:
        model2 = mutate(model,0.05)
        chi2 = chisq(model2,f,xmin)
        if chi2<chi:
            model=model2
            chi=chi2
        pn=old_div(1.0,(1.0+old_div(model['l'],(model['nu']*G))))
        print("#",step,chi,model,model_string(model['nu']*G/pn,model,pn,G,n,fit_params))
#        sys.stdout.flush()
        step+=1

    xs=list(f.keys())
    xs.sort()
    for x in xs:
        y=f.get(x,0)
        print(x,y,fx(x,model))
        z=x
    z=int(z)
    for x in range(z+50000,int(5e6),50000):
        y=f.get(x,0)
        print(x,y,fx(x,model))
        
    pn= old_div(1.0,(1.0+old_div(model['l'],(model['nu']*G))))
    Ne= model['nu']*G/pn
    print("#",pn,Ne,Ne*pn/G)


    if args.outfile:
        f=open(args.outfile,"wt")
        #        f.write(model_string(fit_params))
        if not args.json:
            f.write(model_string(model['nu']*G/pn,model,pn,G,n,fit_params))
            f.write("\n")
        else:
            json.dump(model_json(model['nu']*G/pn,model,pn,G,n,fit_params), f)
        f.close()


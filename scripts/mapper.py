#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import map
from builtins import str
from builtins import range
from builtins import object
import pickle as pickle
debug=False

class node(object):
    def __init__(self,s,x):
        self.s = s
        self.x   = x
#        self.forward_link = False
#        self.reverse_link = False
        self.links = [False,False]
        self.strand = 1
        self.mapped_coord = 0
        self.root = self
        self.leaf = self


    def set_forward_link(self,o):
        self.links[0] = o
        
    def set_reverse_link(self,o):
        self.links[1] = o
        
    def get_forward_link(self):
        return self.links[0]
        
    def get_reverse_link(self):
        return self.links[1] 

    def linkable(self):
        return ((not self.links[0]) or (not self.links[1]))

    def __repr__(self):
        return str((self.s,self.x))#,self.strand,self.mapped_coord))

    def __hash__(self):
        return hash(tuple((self.s,self.x)))

    def __cmp__(self,other):
        if self.__hash__() < other.__hash__():
            return -1
        elif self.__hash__() == other.__hash__():
            return 0
        else:
            return 1

    def __eq__(self,other):
        if self.__hash__() == other.__hash__():
            return True
        else:
            return False

class map_graph(object):
    def __init__(self,scaffold_name_list,scaffold_len_list,minlength=0):
        self.originals_index={} #[]
        self.originals_names={} #list(scaffold_name_list)
        self.debug=False
#        print "names:",self.originals_names
        self.name_to_index = {}
        self.obsoleted={}
        self.roots = {}
        self.length = {}
        self.node={}
        self.alias={}
        self.alias_to_node={}
        self.dont_update_coords=False
        for i in range(len(scaffold_len_list)):
            if scaffold_len_list[i] >= minlength:
                self.originals_names[i]=scaffold_name_list[i]
                self.name_to_index[scaffold_name_list[i]] = i
                ll = [ node(i,1),node(i,scaffold_len_list[i]) ]
                self.length[ (ll[0],ll[1]) ] = scaffold_len_list[i]-1
                self.length[ (ll[1],ll[0]) ] = scaffold_len_list[i]-1
                self.length[ (ll[0].s,ll[0].x) ] = scaffold_len_list[i]
                self.node[ (ll[0].s,ll[0].x) ] = ll[0]
                self.node[ (ll[1].s,ll[1].x) ] = ll[1]
                ll[0].set_forward_link(ll[1])
                ll[1].set_reverse_link(ll[0])
                ll[1].root = ll[0]
                ll[1].leaf = ll[1]
                ll[0].root = ll[0]
                ll[0].leaf = ll[1]
#                ll[0].mapped_coord = 1
                ll[1].mapped_coord = scaffold_len_list[i]-1
                self.roots[ll[0]]=1
                self.originals_index[i]=ll
        #print self.length

    def set_alias(self,node,alias):
        self.alias[node]=alias
        self.alias_to_node[alias]=node
        
    def find_node_by_name(self,name):
        if name in self.alias_to_node:
            print("#alias",name,self.alias_to_node[name])
            return self.alias_to_node[name]
        end=False
        if name[-2:]==".5" or name[-2:]==".3":
            end = name[-1]
            name=name[:-2]
                
        c = name.split("_")
        base = int(c[-1])
        n = "_".join(c[:-1])

        no = self.node[(self.name_to_index[n],base)]
        if end=="3":
            no = no.get_forward_link()
            #print self.print_name(no),self.print_name(ni)


        return no
        
    def find_node_by_pair(self,id,base):
        return self.node[(id,base)]

    def write_map_to_file(self,filename):
        f=open(filename,"wb")
        pickle.dump( self, f)
        f.close()
#        f.write( 
#        f.write(str([ i.serialize() for i in self.originals_index ] ))
#        f.write("\n")
#        f.close()

    def add_oriented_break(self,s,x,strand):
        a,b = self.add_break_in_original(s,x)
        if strand=="+": return a,b
        else: return b,a

    def flip(self,a):
        r,l= a.root,a.leaf
        if r in self.roots : del self.roots[r]
        self.roots[l]=1
        self.reroot(l)

    def add_break_in_original(self,s,x):
        ll = self.originals_index.get(s,False)
        print("add break in",ll,s,x,len(ll))
        if not ll:
            print("warning: tried to add break to unknown scaffold")
            return
        i=len(ll)-1
        while ll[i].x > x:
            print(len(ll),i,x,ll[i])
            i-=1 
        a = node(s,x)
        b = node(s,x+1)
        #        a.root = ll[i].root

        self.node[ (s,x  ) ] = a
        self.node[ (s,x+1) ] = b

        ll.insert(i+1,b)
        ll.insert(i+1,a)

        if ll[i].root in self.roots: del self.roots[ll[i].root]
        a.root = a
        self.roots[a]=1


        if ll[i+3].root in self.roots: del self.roots[ll[i+3].root]
        self.roots[b]=1

        b.root = b
        b.mapped_coord = 0
        a.mapped_coord = 0

#        ll[i+3].root = b

        if not ll[i].get_forward_link()==ll[i+3]:
            print("wtf?",ll[i],ll[i+1],ll[i].get_forward_link(),ll[i+3])
            raise Exception

        ll[i].set_forward_link(a)
        a.set_reverse_link(ll[i])

        self.length[(ll[i],a)]=abs(x-ll[i].x)
        self.length[(a,ll[i])]=abs(x-ll[i].x)

        self.length[(b,a)]=1
        self.length[(a,b)]=1

        b.set_forward_link(ll[i+3])
        ll[i+3].set_reverse_link(b)

        self.length[(ll[i+3],b)]=abs(x+1-ll[i+3].x)
        self.length[(b,ll[i+3])]=abs(x+1-ll[i+3].x)

        self.roots[b]=1
        self.reroot(a.root)
        self.reroot(b.root)

        return(a,b)
#    def remap_coord(self,name,x):
#        n = self.find_node_by_name(name)
#        
#        queue=[n]
#        seen={}
#        while len(queue)>0:
#            s = queue.pop()
#            if s.mapped_coord


    def map_coord_by_name(self,name,x):
        n = self.find_node_by_name(name)
        i=0
        ll = self.originals_index.get(n.s,False) 
        if not ll: 
            print("wtf?", name,x,n)
            exit(0)
        print("##", name,x,n,ll)
        while (i < len(ll)-1) and ll[i+1].x<x:
#            print i, len(ll), ll[i],x,ll[i+1]
            i+=1
        if i%2==1:
            i+=1
#            print i
        xx = ll[i].strand * ( x - ll[i].x ) +  ll[i].mapped_coord + 1

#        print "####",xx,i,ll[i].root,ll[i],ll[i].get_reverse_link(),ll[i].get_forward_link()
#        return((ll[i].s,xx,ll[i].root))
        a,b = ll[i].root.s, ll[i].root.x
        return((a,b),xx)

    def map_coord(self,s,x):
        i=0
        ll = self.originals_index.get(s,False)
        if not ll:
            return((s,1),x)
#        print s,x,ll
        while (i < len(ll)-1) and ll[i+1].x-1<=x:
#            print i, len(ll), ll[i],x,ll[i+1]
            i+=1
        if self.debug: print(s,x,i)

#XXX why was this ever done?
#        if i%2==1:
#            if self.debug: print "XXXX"
#            i+=1
#            print i
        if self.debug: print(i, len(ll))
        xx = ll[i].strand * ( x - (ll[i].x -1) ) +  ll[i].mapped_coord 

        if xx<0 : 
            print("wtf?", ll[i].strand * ( x - (ll[i].x-1) ) +  ll[i].mapped_coord, "=", ll[i].strand ,"* (", x," -", (ll[i].x-1)," ) + ", ll[i].mapped_coord, "ll[i].strand * ( x - (ll[i].x-1) ) +  ll[i].mapped_coord")
            print(self.print_name(ll[i]))
            print("x=",x)
            print("i",i)
            print("s=",s)
            print("ll[i]=",ll[i])
            print("ll[i].x=",ll[i].x)
            print("ll=",ll)
            print("ll[i].mapped_coord=",ll[i].mapped_coord)
            #raise Exception 


#        print "####",xx,i,ll[i].root,ll[i],ll[i].get_reverse_link(),ll[i].get_forward_link()
#        return((ll[i].s,xx,ll[i].root))
        a,b = ll[i].root.s, ll[i].root.x
        return((a,b),xx)

    def find_original_coord(self,node,y,gapslop=0):
        #ll = self.originals_index(node.s,False)
        #if not ll: print "wtf?",s,x,"find_originals_coord"
        onode=node
        start = node.mapped_coord
        queue = [node]
        last=False
        visited={}
        b1=False
        b2=False
        ox=-1
    
        while len(queue)>0:
            n = queue.pop(0)
            if not visited.get(n,False):
                if last:
                    print("visit:",last,n,last.mapped_coord,n.mapped_coord,y)
                    if ( last.mapped_coord <= y  and y < n.mapped_coord):
                        if (last.s == n.s):  

                            if (y-last.mapped_coord)<gapslop:
                                print("bump back:",y,last2,last)
                                return(last2,last,"gap")
                            elif (n.mapped_coord - y)<gapslop:
                                print("bump ahead:",y,n.mapped_coord+1)
                                y=n.mapped_coord+1
                            else:
                                strand="+"
                                if last.x>n.x : strand="-"
                                s = last.s
                                if strand=="+":
                                    ox=last.x+y-last.mapped_coord
                                else:
                                    ox=last.x-(y-last.mapped_coord)
                                print("forward",last, n, last.mapped_coord,"<=",y,"<", n.mapped_coord, self.print_name(last), self.print_name(n),ox)
    #                            return(self.originals_names[last.s],ox, strand )
                                return(last.s,ox, strand )
                        else:
    #                       print "gap forward", n, last, n.mapped_coord, last.mapped_coord, self.print_name(n), self.print_name(last)
                            print("gap forward", last, n, last.mapped_coord,"<=",y,"<", n.mapped_coord, self.print_name(last), self.print_name(n))
    #                        print "gap",last,n
                            return(last,n,"gap")
#                            return((last.s,last.x),(n.s,n.x),"gap"+strand)

                visited[n]=True
                if n.get_forward_link(): queue.append(n.get_forward_link())
                if n.get_reverse_link(): queue.append(n.get_reverse_link())
#                if n.get_reverse_link() == last: 
#                    #n.strand = 1
#                    pass
#                else:
#                    pass
                    #n.strand = -1
                last2 = last
                last  = n

#        if b1 and b2:
#            print node,node.mapped_coord ,node.leaf,node.leaf.mapped_coord, y,b1,b2,b1.mapped_coord,b2.mapped_coord, b1.s,b1.x, b2.s,b2.x
#        else:
#            print node,node.mapped_coord ,node.leaf,node.leaf.mapped_coord, y,b1,b2 #,b1.mapped_coord,b2.mapped_coord
        print("mapper wtf?",self.print_name(onode),self.originals_names[onode.s],y,self.print_name(n),n.mapped_coord,ox,visited)
        return(n,ox)

#(137628, 1) (137628, 3088) 3087 3088 False False



    def originals(self,a):
        oris = {}
        queue=[a]
        visited={}
        mapped_coord = 0
        a.mapped_coord = 0
        visited[a]=True
        if a.get_forward_link(): 
            queue.append(a.get_forward_link())
            a.strand = 1
        if a.get_reverse_link(): 
            queue.append(a.get_reverse_link())
            a.strand = -1
        last = a
        root = a
        a.root = a
        while 0<len(queue):
            n = queue.pop(0)
            oris[n.s]=1
            if not visited.get(n,False):
                visited[n]=True            
                if n.get_forward_link(): queue.append(n.get_forward_link())
                if n.get_reverse_link(): queue.append(n.get_reverse_link())
                if n.get_reverse_link() == last: 
                    n.strand = 1
                else:
                    n.strand = -1
                last = n
        return(list(oris.keys()));
 
    def print_name(self,n):
#        print self.originals_names
#        print n,self.originals_names[n[0]],"_",str(n[1])

        if n in self.alias: return self.alias[n]

        if type(n)==type((1,1)):
            return self.originals_names[n[0]]+"_"+str(n[1])
        else:
            return self.originals_names[n.s]+"_"+str(n.x)

    def nodepair_list(self,n):
        lasti=-1
        lastx=-1
        lastn=False
        for nn in self.nodes_dfs(n):
            if lasti==nn.s:
                yield( ( lastn , nn )  )
            lastn=nn
            lasti=nn.s

    def region_list(self,n):
        lasti=0
        lastx=-1
        for nn in self.nodes_dfs(n):
#            print nn
            if lasti==nn.s:
                if lastx<nn.x:
                    strand=1
                else:
                    strand=-1
                range_start = min(lastx,nn.x)
                range_end = max(lastx,nn.x)
                if range_start>0 and range_end>0:
                    yield (self.originals_names[nn.s]+":"+str(range_start)+"-"+str(range_end),strand,self.originals_names[nn.s],range_start,range_end,min(lastm,nn.mapped_coord),max(lastm,nn.mapped_coord))
            lasti=nn.s
            lastx=nn.x
            lastm=nn.mapped_coord

    def nodes_dfs(self,n,colors={}):
#        yield(n)
        colors={}
        colors[n]=1
        queue = [n]
        while len(queue)>0:
            n = queue.pop(0)
#            print "x#",self.print_name(n)
            yield(n)
            c1 = n.get_forward_link()
 #           print c1
            if c1 and c1 not in colors:
                queue.append(c1)
                colors[c1]=1
    #            self.nodes_dfs(c1,colors)

            c2 = n.get_reverse_link()
  #          print c2
            if c2 and c2 not in colors:
                colors[c2]=1
                queue.append(c2)
    #            self.nodes_dfs(c2,colors)

    def write_dot(self,filename):
        f=open(filename,"wt")
        for r in list(self.roots.keys()):
            f.write(" -- ".join( map(str,self.nodes_dfs(r)) ))
#            for n in self.nodes_dfs(r):
#                f.write( str(n))
#                f.write(" -- ")
#            f.write(str(r))
            f.write("\n")
        f.close()

    def reroot(self,a):
        debug=False
        if debug: print("reroot")
        total_len = 0
        queue=[a]
        visited={}

        if not a == a.root:
            if a.root in self.roots:
                del self.roots[a.root]
                self.roots[a]=1
    
        mapped_coord = 0
        a.mapped_coord = 0
        visited[a]=True
        if a.get_forward_link(): 
            queue.append(a.get_forward_link())
            a.strand = 1
        if a.get_reverse_link(): 
            queue.append(a.get_reverse_link())
            a.strand = -1
        last = a
        root = a
        a.root = a
        if debug: print("root",a)
        if debug: print(">>>>",a)
        while 0<len(queue):
#            print queue,visited
            n = queue.pop(0)
            if not visited.get(n,False):
                visited[n]=True            
                n.root = root
                if n.get_forward_link(): queue.append(n.get_forward_link())
                if n.get_reverse_link(): queue.append(n.get_reverse_link())
                if n.get_reverse_link() == last: 
                    n.strand = 1
                else:
                    n.strand = -1
                mapped_coord += self.length[(last,n)]
                n.mapped_coord = mapped_coord
                if debug: print(">>>>",n,self.length[(last,n)],mapped_coord,root)
                total_len += self.length[(last,n)]
                last = n
        for k in list(visited.keys()):
            k.leaf = last
        self.length[a] = total_len

        if debug: print("done rerooting")

    def update_roots(self):
#        import set
        rr=set(self.roots.keys())
        ob=set(self.obsoleted.keys())
        self.roots={}
        for rrr in rr.difference(ob):
            self.roots[rrr]=1
#        self.roots = list(rr.difference(ob))
#        self.roots = [ r for r in self.roots if not self.obsoleted.has_key(r) ]
        self.obsoleted = {}

    def unjoin(self,a,b):
        if   a.get_forward_link()==b: 
            a.set_forward_link(False)
            #print "unlink a forward"
        elif a.get_reverse_link()==b: 
            a.set_reverse_link(False)
            #print "unlink a reverse"
    
        if   b.get_forward_link()==a: 
            b.set_forward_link(False)
            #print "unlink b forward"
            
        elif b.get_reverse_link()==a: 
            b.set_reverse_link(False)
            #print "unlink b reverse"

        self.reroot(a)
        self.reroot(b)

    def add_join(self,a,b,g=1):
        debug=False
        if debug: print("add join",a,b,g)

        self.length[(a,b)] = g+1
        self.length[(b,a)] = g+1
        if debug: print(a,b,a.root,b.root)

        if a.root == b.root:
            print("Tried to join two nodes already having the same root.", a,b,a.root,b.root,a.leaf,b.leaf)
            return

        newroot = a.root
        if newroot == a:
#            print "##root clash"
            newroot = b.root
        if newroot == b:
            newroot = a.leaf

        if debug: print("length",a,b,g)
        if debug: print("length",b,a,g)

        if not newroot == a.root:
#            self.obsoleted[a.root]=True
            if a.root in self.roots: del self.roots[a.root]
        if not newroot == b.root:
#            self.obsoleted[b.root]=True
#            del self.roots[b.root]
            if b.root in self.roots: del self.roots[b.root]

        #  --a--> --b-->
        if (not a.get_forward_link()) and (not b.get_reverse_link()):
            if debug: print(" --a--> --b--> ")
            a.set_forward_link(b)
            b.set_reverse_link(a)
  #          self.roots.append(a.root)
   #         self.reroot(a.root)

        #  --a--> <--b--
        elif (not a.get_forward_link()) and (not b.get_forward_link()):
            if debug: print(" --a--> <--b-- ")
            a.set_forward_link(b)
            b.set_forward_link(a)
 #           self.roots.append(a.root)
  #          self.reroot(a.root)

        #  <--a-- --b-->
        elif (not a.get_reverse_link()) and (not b.get_reverse_link()):
            if debug: print(" <--a-- --b--> ")
            a.set_reverse_link(b)
            b.set_reverse_link(a)
#            self.roots.append(a.leaf)
 #           self.reroot(a.leaf)
            
        #  <--a-- <--b--
        elif (not a.get_reverse_link()) and (not b.get_forward_link()):
            if debug: print(" <--a-- <--b-- ")
            a.set_reverse_link(b)
            b.set_forward_link(a)
#            self.roots
#            self.roots.append(b.root)
#            self.reroot(b.root)
        else:
            print("Tried to join to edges that are already joined",a,b)
            exit(0)

        self.roots[newroot]=1
        if not self.dont_update_coords:
            print("#rerooting")
            self.reroot(newroot)
        debug=False


if __name__=="__main__":

    m = map_graph( ["a","b","c"],[10000,10000,10000] )
    m.add_break_in_original(1,5000)

    for r in list(m.roots.keys()):
        print("scaffold:",r, "--" ,r.leaf)


    print("#",m.map_coord(1,4001),4001)
    print("#",m.map_coord(1,6001),1001)

    m.add_join( m.originals_index[1][1],m.originals_index[1][2] )
    #m.add_join( m.originals_index[0][1],m.originals_index[2][0] )

    for r in list(m.roots.keys()):
        print("scaffold:",r, "--" ,r.leaf)

    print("#",m.map_coord(1,4001),4001)
    print("#",m.map_coord(1,6001),6001)
    print("#",m.map_coord(0,1001),1001)
    print("#",m.map_coord(2,1001),1001)

    m.add_join( m.originals_index[1][3],m.originals_index[2][0] )

    for r in list(m.roots.keys()):
        print("scaffold:",r, "--" ,r.leaf)

    print("#",m.map_coord(1,4001),4001)
    print("#",m.map_coord(1,6001),6001)
    print("#",m.map_coord(0,1001),1001)
    print("#",m.map_coord(2,1001),11001)

    m.write_map_to_file("test_map.txt")


    m.add_join( m.originals_index[1][0],m.originals_index[0][0] )

    for r in list(m.roots.keys()):
        print("scaffold:",r, "--" ,r.leaf)

    print("#",m.map_coord(1,4001),4001)
    print("#",m.map_coord(1,6001),6001)
    print("#",m.map_coord(0,100) ,29900)
    print("#",m.map_coord(2,1001),11001)

    f = open("test_map.txt")
    m2 = pickle.load(f)
    print(m2)

    for r in list(m2.roots.keys()):
        print("scaffold:",r, "--" ,r.leaf)




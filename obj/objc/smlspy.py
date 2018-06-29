#!/usr/bin/env python3
# coding=utf-8

import string
import os
import re
import numpy as np
# import scipy as sp
# import matplotlib as mp

class smlspy:
    """
    the module is for Smiles, containing the conversion a Smiles format to
    a relation matrix, and additional function based on this.

    the functions as follow:
    
    Class smlspy:
        __init__(): to get a smlspy object
        sml():
        smat():
        dous():
        cmpr():
    """
    am={'h':1.001,'c':12.01,'n':14.01,'o':16.00,'f':19.00,'p':30.97,'s':32.07}
    od={'c':1,'n':2,'o':3,'f':4,'p':5,'s':6,'x',7}
    hy1={'h':1,'c':4,'n':3,'o':2,'f':1,'p':3,'s':2}
    hy2={'h':1,'c':4,'n':5,'o':2,'f':1,'p':5,'s':6} # highest valence
    
    errinfo = 0; # show the level of error

    def __init__(self, expr):
        self.state=0
        self.expr=expr
        self.N=len(expr)
        self.sml='' # sml format sequence
        self.std='' # standard format sequence
        self.stdmap=[] # get the std sequence
        self.map=[] # record the index number of sml in initial expr, type int list
        self.amap=[] # reversed map can be useful
        self.clip=[] # record the clip
        self.frag=[] # record the position of each fraction

        # the major of two value map & sml, will be calculate at first
        cnt=0
        for i in range(len(expr)):
            if expr[i].isalpha():
                self.map.append(i)
                self.amap.append(cnt)
                cnt+=1
            else:
                self.amap.append(-1)
        self.n=len(self.map)

        for i in self.map:
            self.sml+=self.expr[i]

    def getclip(self):
        cnt=1
        ptn=r'[A-Za-z]'+str(cnt)+'[A-Za-z\(\)=][A-Za-z0-9\(\)=]*[A-Za-z]'+str(cnt) # maybe err
        
        srch=re.search(ptn,self.expr)
        while srch:
            (b,e)=srch.span()
            b=int(b);e=int(e)-1
            
            while not self.expr[e].isalpha(): # find the atom position of the end
                e-=1
            
            self.clip.append([b, e]) # add the clipping head and tail to the clip list
            cnt+=1
            ptn=r'[A-Za-z]'+str(cnt)+'[A-Za-z\(\)=][A-Za-z0-9\(\)=]*[A-Za-z]'+str(cnt)
            srch=re.search(ptn,self.expr)
        return len(self.clip)
    
    def getfrag(self):
        exprx=self.expr
        ptn=r'\(=?[A-Za-z]' # att at C(=CC) structure!
        srch=re.search(ptn, exprx)
        while srch:
            (b,e)=srch.span()
            b=int(b)
            ip=b # use for search bracket
            exprx=exprx[0:b]+'{'+exprx[b+1:-1] # change state of exprx

            if exprx[b+1]=='=': # b point at the (atom
                b+=2
            else:
                b+=1
            
            next = True
            bra=0
            while next: # search pre-atom
                if exprx[ip]=='(' or exprx[ip]=='{':
                    bra+=1
                if exprx[ip]==')' or exprx[ip]=='}':
                    bra-=1
                
                if bool(bra>=1) & bool(re.search(r'[A-Za-z][0-9]*$',exprx[0:ip])):
                    ip-=1
                    while exprx[ip].isdigit():
                        ip-=1
                    self.frag.append([ip,b])
                    next=False
                ip-=1
                if ip<=0:
                    self.errinfo=1
            ptn=r'\(=?[A-Za-z]' # att at C(=CC) structure
            srch=re.search(ptn, exprx)
        
        ptn=r'\)=?[A-Za-z]' # att at C(C)=CC structure
        srch=re.search(ptn, exprx)
        while srch:
            (b,e)=srch.span()
            b=int(b)
            ip=b # use for search bracket
            exprx=exprx[0:b]+'}'+exprx[b+1:-1] # change state of exprx

            if exprx[b+1]=='=': # b point at the (atom
                b+=2
            else:
                b+=1
            
            next = True
            bra=0
            while next: # search pre-atom
                if exprx[ip]=='(' or exprx[ip]=='{':
                    bra+=1
                if exprx[ip]==')' or exprx[ip]=='}':
                    bra-=1
                
                if bool(bra>=0) & bool(re.search(r'[A-Za-z][0-9]*$',exprx[0:ip])):
                    ip-=1
                    while exprx[ip].isdigit():
                        ip-=1
                    self.frag.append([ip,b])
                    next=False
                ip-=1
                if ip <0:
                    self.errinfo=1
                    next = False
            ptn=r'\)=?[A-Za-z]' # att at C(=CC) structure!
            srch=re.search(ptn, exprx)


    def ckbond(self,i,j):
        if max(i,j)>=self.n or min(i,j) < 0:
            errinfo=-1 # index over the boundary
            return 0;
        if i==j: # to judge the hyd-cls of the atom
            spn=self.hy1[self.expr[self.map[i]].lower()]
            if self.expr[self.map[i]].islower():
                spn-=1
            if self.expr[self.map[i-1]]=='=':
                spn-=1
            if self.expr[self.map[i+1]]=='=':
                spn-=1
            for tm in self.frag:
                if tm[0]==self.map[i]:
                    if self.expr[tm[1]-1]=='=':
                        spn=-1
                    if self.expr[tm[1]-2]=='=':
                        spn-=1
            if spn<=0:
                spn=spn+self.hy2[self.expr[self.map[i]].lower()]-self.hy1[self.expr[self.map[i]].lower()]
            if spn>0:
                return spn
            else:
                self.errinfo=-3
                self.waring() # warning------------------------------------3
                return 0;
        if i > j: # get the type of bond
            k=i;i=j;j=k;
        if self.map[i]+1==self.map[j]:
            if self.expr[self.map[i]].islower() & self.expr[self.map[j]].islower():
                return 2;
            else:
                return 1;
        elif self.map[i]+2==self.map[j]:
            if self.expr[self.map[i]+1]=='=':
                if self.expr[self.map[i]].islower() & self.expr[self.map[j]].islower():
                    return 3;
                else:
                    return 2;
        elif self.map[i]+3==self.map[j]:
            if bool(self.expr[self.map[i]+1]=='=') & bool(self.expr[self.map[i]+2]=='='):
                return 3;
        else:
            return 0;

    def lth(self): # not useful
        lthmt=np.mat(np.zero(self.n))
        return 0

    def wht(self): # not useful
        wht=[0]*self.n
        for i in range(self(n)):
            for j in range(i+1,self.n):
                lth= self.lth(i,j)
                wht[i]+=(10*0.1**lth)*self.am[self.expr[self.map[j]].lower()]
                wht[j]+=(10*0.1**lth)*self.am[self.expr[self.map[i]].lower()]
        return 0;

    def sym(self):


    def getstdmap(self):
        smlmt=self.smlmat()
        seq=[0]*self.n
        for i in range(self.n):
            sumi=self.am[self.expr[self.map[i]].lower()]*100
            for j in range(i+1,self.n):
                sumi+=smlmt[i,j]
            seq[i]=sumi
        for i in range(self.n):
            idx=seq.index(max(seq))
            self.stdmap.append(idx)
            seq[idx]=0
        return 0


    def smlmat(self):
        if self.clip==[]:
            self.getclip()
        if self.frag==[]:
            self.getfrag()
        smlmt=np.mat(np.eye(self.n))
        for i in range(self.n):
            for j in range(i+1,min(i+3,self.n)):
                if self.expr[self.map[j]]=='(' or self.expr[self.map[j]]==')':
                    break
                if  self.ckbond(i,j)==1:
                    smlmt[i,j]=1
                elif self.ckbond(i,j)==2:
                    smlmt[i,j]=2
                else:
                    smlmt[i,j]=0
        for tm in self.clip:
            smlmt[self.amap[tm[0]],self.amap[tm[1]]]=1
        for tm in self.frag:
            smlmt[self.amap[tm[0]],self.amap[tm[1]]]=1 # or 2?
        return smlmt

    def dous(self):
        if self.clip==[]:
            self.getclip()
        d=0
        for i in range(self.N):
            if self.expr[i].islower():
                d+=1
            if self.expr[i]=='=':
                if self.expr[i-1]=='=' or self.expr[i+1]=='=':
                    d+=1.5
                else:
                    d+=2
        d=int(len(self.clip)+d/2)
        return d

    def isequ(self, ob):
        if self.n != ob.n:
            return False
        smlmt1=self.smlmat()
        self.getstdmap()
        smlmt2=ob.smlmat()
        ob.getstdmap()

        for i in range(self.n):
            si=self.stdmap[i]
            obi=ob.stdmap[i]
            if self.expr[self.map[si]].lower() != ob.expr[ob.map[obi]].lower():
                return False
            for j in range(i,self.n):
                sj=self.stdmap[j]
                obj=ob.stdmap[j]
                if self.ckbond(si,sj) != ob.ckbond(obi,obj):
                    return False
        return True

if __name__ == '__main__':
    test=smlspy('C1CC((C=C2)(OC))C2CCC1')
    test2=smlspy('C1CC((cc2)OC)C2CCC1')
    print(type(test))
    print('the original seq is:\n\t',test.expr,'\n\t',test2.expr)
    print('origin num is: ',test.N,test2.N)
    print('active num is: ',test.n,test2.n)
    
    print('the map is:\n ',test.map,test2.map)
    print('the amap is:\n ',test.amap,test2.amap)
    print('the sml seq is:\n\t',str(test.sml),str(test2.sml))
          
    print(test.frag)
    print(test.smlmat())
    print(test2.smlmat())
    print(test.isequ(test2))

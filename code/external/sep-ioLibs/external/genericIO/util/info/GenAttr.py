#!/usr/bin/env python3
import argparse
import genericIO
import GenJob
import GenSplit
import numpy as np
import numba
import copy
import math
import threading


import copy

class helix2cart:
    """Convert to and from cartesian and helix space"""
    def __init__(self,nd):
        """Initialize conversion

             nd - Data dimensions"""
        if not isinstance(nd,list):
            raise Exception("Expecting nd to be a list")

        self._ndim=copy.deepcopy(nd)
        self._b=[1]
        sz=1
        for n in self._ndim:
            if not isinstance(n,int):
                raise Exception("Expecting a list of ints")
            sz=sz*n
            self._b.append(sz)


    def toCart(self,hlx):
        cart=[0]*len(self._ndim)
        lft=hlx
        for i in range(len(self._ndim)-1,-1,-1):
            cart[i]=int(lft/self._b[i])
            lft-=cart[i]*self._b[i]
        return cart


    def toHelix(self,cart):
        """Convert from cartesian space to helix space"""
        if len(cart) != len(self._ndim):
            raise Exception("Expecting cart to be same size as data")
        hlx=0
        for i in range(len(self._ndim)):
            if not isinstance(cart[i],int):
                raise Exception("Expecting cart to be a list of ints")
            hlx+=cart[i]*self._b[i]
        return hlx

class attrJob(GenJob.regSpace):
    def __init__(self,inputType):
        """Intialize object

            inputType - Input type
        """
        super().__init__(self.calcStats,0,0,inputType=inputType,hasOutput=False)
        self._lock=threading.Lock() 
        self._mn=1e99
        self._mx=-1e99
        self._imin=-1
        self._imax=-1
        self._sm=0
        self._sqs=0
        self._nzero=0
    def calcStats(self,ina,dummy):
        """Convert a buffer from one type to another

        ina - Input vector
        dummy - Dummy argument
        """

        n123=ina.getHyper().getN123()
        inN=np.reshape(ina.getNdArray(),(n123,))

        if self._inputType=="dataComplex" or self._inputType=="dataComplexDouble":
            mn,imin,mx,imax,sm,sqs,nzero=calcComplexStats(inN)
        else:
           mn,imin,mx,imax,sm,sqs,nzero= calcRealStats(inN)
        self._lock.acquire()
        if mn< self._mn:
            self._mn=mn
            self._imin=imin
        if mx > self._mx:
            self._mx=mx
            self._imax=imax 
        self._sm+=sm
        self._sqs+=sqs
        self._nzero+=nzero
        self._lock.release()

        
@numba.jit(nopython=True, parallel=True,locals={'sm': numba.float64,"sqs":numba.float64,"nzero":numba.int64})
def calcRealStats(inA):
    """
      Return min,max,sum,sumsq,nzeros
      """
    nzero=0
    sm=inA[0]
    sqs=inA[0]
    mn=inA[0]
    mx=inA[0]
    imin=0
    imax=0
    if inA[0]!=0: 
        nzero+=1
    for i in range(1,inA.shape[0]):
        if inA[i]  < mn:
            mn=inA[i]
            imin=i
        if inA[i] >mx:
            mx=inA[i]
            imax=i
        sm+=inA[i]
        sqs+=inA[i]*inA[i]
        if inA[i] !=0:
            nzero+=1
    return mn,imin,mx,imax,sm,sqs,nzero

@numba.jit(nopython=True,locals={'amp':numba.float64,'sm': numba.float64,"sqs":numba.float64,"nzero":numba.int64})
def calcComplexStats(inA):
    """
      Return min,max,sum,sumsq,nzeros
      """
    nzero=0
    amplitude=np.absolute(inA[0])
    sm=amplitude
    sqs=amplitude
    mn=amplitude
    mx=amplitude
    imin=0
    imax=0
    if amplitude !=0: 
        nzero+=1
    for i in range(1,inA.shape[0]):
        amplitude=numpy.absolutde(inA[0])
        if amplitude  < mn:
            mn=amplitude
            imin=i
        if amplitude >mx:
            mx=amplitude
            imax=i
        sm+=amplitude
        sqs+=amplitude*amplitude
        if amplitude !=0:
            nzero+=1
        
    return mn,imin,mx,imax,sm,sqs,nzero
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Print data attributes')
    parser.add_argument('input', metavar='Input', type=str,
                        help='Input file')                 
    parser.add_argument("--io", type=str,choices=[@GEN_IO_TYPES@], help='IO type. Defaults to defaultIO')
    parser.add_argument("--memory",type=int,help="Memory in terms of GB",default=20)
    parser.add_argument("--want",type=str,choices=["all","min","max","mean","rms","norm","short"],
    default="all")
    parser.add_argument("--print_pct",type=float,help="Print progress every X pct (above 100 means no printing)",default=101)
    args = parser.parse_args()

    ioIn=genericIO.defaultIO

    if args.io:
        ioIn=genericIO.io(args.io)
    


    inFile=ioIn.getRegFile(args.input)
    job=attrJob(inFile.getStorageType())
    job.setCompleteHyperOut(inFile.getHyper())
    job.setInputFile(inFile)
    split=GenSplit.serialRegSpace(job, args.memory)
    split.loop(args.print_pct)
    
    hx=helix2cart(inFile.getHyper().getNs())

    n123=inFile.getHyper().getN123()
    if args.want=="all":
        print("rms=%f"%(math.sqrt(job._sqs/n123)))
        print("mean=%f"%(job._sm/n123))
        print("norm=%f"%(math.sqrt(job._sqs)/n123)
        print("maximum value =%f at "%job._mx,hx.toCart(job._imax))
        print("minimum value =%f at "%job._mn,hx.toCart(job._imin))
        print("number of nonzero sammples=",job._nzero)
        print("total number of samples=",n123)
    elif args.want=="min":
        print(job._min)
    elif args.want=="max":
        print(job._max)
    elif args.want=="mean":
        print(job._sm/n123)
    elif args.want=="rms":
        print(math.sqrt(job._sqs/n123))
    elif args.want=="norm":
        print(math.sqrt(job._sqs)/n123)
    else:
        print(job._nzero/n123)









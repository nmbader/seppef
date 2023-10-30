#!/usr/bin/env python3
import argparse
import genericIO
import GenJob
import GenSplit
import Hypercube
import math
import numpy as np
from numba 

class sincTable:
    def __init__(self,ntab,nlen):
        if nlen != int(nlen/2)*2:
            raise Exception("nlen must be divisible by two")
        self._dtab=1./float(ntab)
        self._tab=np.zeros((ntab,nlen))
        b=np.zeros((nlen,))
        c=np.zeros((nlen,))
        work=np.zeros((nlen,))
        mkTable(self._tab,b,c,work)
    
    def getSinc(self,val):
        """Return a numpy array that best represents val. val has to be between 0 and 1"""
        if val <0. or val >=1.:
            raise Exception("val must greater or equal to zero and less than one")
        ival=int(round(val/self._dtab))
        return self._tab[ival,:]

@numba.jit(nopython=True,locals={"c": numba.float64, "e": numba.float64, "v": numba.float64, "w": numba.float64, "bot": numba.float64})
def toep(r,f,g,a):
    a[0]=1
    v=r[0]
    f[0]=g[0]/r[0]
    for j in range(1,r.shape[0]):
        e = 0
        a[j] = 0
        f[j] = 0.
        for i in range(r.shape[0]):
            e+=a[i]*r[j-i]
        c= e/v
        v-= e*c 
        jh = int(j/2)
        for i in range(jh):
            bot= a[j-i] - c* a[i]
            a[i] -= c*a[j-i]
            a[j-i] =bot
        w=0
        for i in range(r.shape[0]):
            w+= f[i] * r[j-i]
        c= (g[j]-w)/v
        for i in range(j+1):
            f[i]+= c* a[j-i]
    return f



@jit(nopython=True,parallel=True)
def mkTable(tab,b,c,work):
    dtab=1./float(tab.shape[0])
    pi = 3.141592654;
    pi2 = pi*2;
    snyq = 0.5;
    snat = snyq*(0.066+0.265*math.log(double(tab.shape[1]))
    s0 = 0.0;
    ds = (snat-s0)/(tab.shape[1]*2-1);
    for itab in numba.prange(tab.shape[0]):
        d=itab*dtab
        eta = tab.shape[1]/2-1.0+d;
        for j in tab.shape[1]:
            s=s0
            b[j]=0.
            c[j]=0.
            while s < snat:
                b[j]+=math.cos(pi2*s*j)
                c[j]+=math.cos(pi2*s*(eta-j))
                s+=ds
        tab[itab,:]=toep(b,sinc,c,work)
            

class interpJob(GenJob.regSpace):
    def __init__(self,inFile,outFile,interpType=2,nsincLen=10):
        """Intialize object
            inFile - Input file
            outFile - Output file
            interpType - 0 nearest neighbor, 1 linear, 2 sinc interpolation
            nsincLen - Length of sinc. If sinc interpolation is used
      
        """
        if interpType==2:
            self._sincTable=sincTable(10000,10)
        
        if self._interpType=interpType

        if inFile.storage != outFile.storage:
            raise Exception("Storage of input and output not the same")
        
        super().__init__(self.interpOp,0,0,inputType=inFile.storage ,outputType=outFile.storage)
        self._events=events
    
    def interpOp(self,ina,outa):
        """Convert a buffer from one type to another

        ina - Input vector
        outa - Output vector
        """
        nin=np.asarray([1]*6),dtype=np.int32)
        nout=np.asarray([1]*6),dtype=np.int32)
        ooin=np.asarray([0.]*6),dtype=np.float32)
        oout=np.asarray([0.]*6),dtype=np.float32)
        din=np.asarray([1.]*6),dtype=np.float32)
        dout=np.asarray([1.]*6),dtype=np.float32)
        axesIn=ina.getHyper().axes
        axesOut=outa.getHyper().axes
        for i in range(len(axesIn)):
            nin[i]=axesIn[i].n 
            oin[i]=axesIn[i].o
            din[i]=axesIn[i].d
        for i in range(len(axesOut)):
            nout[i]=axesOut[i].n 
            oout[i]=axesOut[i].o
            dout[i]=axesOut[i].d

        if self._interpType==0:
            self.nearestInterp(ina,nin,oin,din,outa,nout,oout,dout)
        elif self._interpType==1:
            self.linearInterp(ina,nin,oin,din,outa,nout,oout,dout)
        else:
            self.sincInterp(ina,nin,oin,din,outa,nout,oout,dout)


    def nearestInterp(ina,nin,oin,din,outa,nout,oout,dout):
        """Do nearest neighbor interpolation"""
        for i6 in range(nout[5])
            out6=oout[5]+dout[5]*i6
            j6=max(0,min(nin[5]-1,int(round((out6-oout[5])/dout[5]))))
            for i5 in range(nout[4]):
                out5=oout[4]+dout[4]*i5
                j5=max(0,min(nin[4]-1,int(round(out5-oout[5]/dout[4]))))  
                for i4 in range(nout[3])
                    out4=oout[3]+dout[3]*i4
                    j4=max(0,min(nin[3]-1,int(round((out4-oout[3])/dout[3]))))
                    for i3 in range(nout[2]):
                        out5=oout[2]+dout[2]*i3
                        j3=max(0,min(nin[2]-1,int(round(out3-oout[2]/dout[2]))))            
                        for i2 in range(nout[1])
                            out2=oout[1]+dout[1]*i2
                            j2=max(0,min(nin[1]-1,int(round((out2-oout[1])/dout[1]))))
                            for i5 in range(nout[4]):
                                out1=oout[0]+dout[0]*i1
                                j1=max(0,min(nin[0]-1,int(round(out1-oout[0]/dout[0]))))
                                outa[i6,i5,i4,i3,i2,i1]=ina[j6,j5,j4,j3,j2,j1] 
                                
                                 
def fill(ar,mag,f1,f2,f3,f4,f5,f6,k1,k2,k3,k4,k5,k6):
    for i6 in k6:
        j6=i6-f6
        if j6 >=0 and j6 < ar.shape[0]:
            for i5 in k5:
                j5=i5-f5
                if j5 >=0 and j5 <ar.shape[1]:
                    for i4 in k4:
                        j4=i4-f4
                        if j4 >=0 and j4 < ar.shape[2]:
                            for i3 in k3:
                                j3=i3-f3
                                if j3>=0 and j3 < ar.shape[3]:
                                    for i2 in k2:
                                        j2=i2-f2
                                        if j2>=0 and j2<ar.shape[4]:
                                            for i1 in k1:
                                                j1=i1-f1
                                                if j1 >=0 and j1<ar.shape[5]:
                                                    ar[j6,j5,j4,j3,j2,j1]=mag
    

def hyperFromArgs(args):
    """Create a hypercube from arguments"""
    axs=[]
    axs.append(Hypercube.axis(n=args.n1,o=args.o1,d=args.d1,label=args.label1))
    axs.append(Hypercube.axis(n=args.n2,o=args.o2,d=args.d2,label=args.label2))
    axs.append(Hypercube.axis(n=args.n3,o=args.o3,d=args.d3,label=args.label3))
    axs.append(Hypercube.axis(n=args.n4,o=args.o4,d=args.d4,label=args.label4))
    axs.append(Hypercube.axis(n=args.n5,o=args.o5,d=args.d5,label=args.label5))
    axs.append(Hypercube.axis(n=args.n6,o=args.o6,d=args.d6,label=args.label6))

    return Hypercube.hypercube(axes=axs)

class event:
    """Create an event for the spike"""
    def __init__(self,ievent,hyper,storage,k1,k2,k3,k4,k5,k6,mag):
        """Event 
            hyper - Hypercube of output
            storage - Storage type
            ievent - Event number
            k1   - Location first axis (fortran)
            k2   - Location second axis (fortran)
            k3   - Location third axis (fortran)
            k4   - Location fourth axis (fortran)
            k5   - Location fifth axis (fortran)
            k6   - Location sixth axis (fortran)
  
            mag  - Magnitude for event"""
        if ievent>= len(mag):
            raise Exception("Event number larger than the magnitude list")


        if storage=="dataInt" or storage=="dataInt" or storage=="dataShort":
            try:
                self._mag=int(mag[ievent])
            except:
                raise Exception("Trouble converting to int:",self._mag[ievent])
            if storage=="dataByte":
                if self._mag < 0 or self._mag >255:
                    raise Exception("For byte expecting values between 0 and 255")
            if storage=="dataShort":
                if self._mag < -32767 or self._mag >32767:
                    raise Exception("For short expecting values between -32767 and 32767")
            if storage=="dataInt":
                if self._mag < -2147483648 or self._mag >2147483648:
                    raise Exception("For short expecting values between -2147483648 and 2147483648")
        elif storage=="dataFloat" or storage=="dataDouble":
            try:
                self._mag=float(mag[ievent])
            except:
                raise Exception("Trouble converting to float:",self._mag[ievent])
        elif storage=="dataComplex" or storage=="dataComplexDouble":
            try:
                self._mag=complex(mag[ievent])
            except:
                raise Exception("Trouble converting to complex:",self._mag[ievent])   
        else:
            raise Exception("Unknown data type ",storage)             

        if ievent >= len(k1):
            self._k1=np.asarray(range(hyper.axes[0].n),dtype=np.int32)
        else:
            k1[ievent]=int(k1[ievent])
            if k1[ievent]<=0 or k1[ievent] > hyper.axes[0].n:
                raise Exception("Illegal value of k1=%d must be between 1 and %d"%(k1[ievent],hyper.axes[0].n))
            self._k1=np.asarray([k1[ievent]-1],dtype=np.int32)

        if ievent >= len(k2):
            self._k2=np.asarray(range(hyper.axes[1].n),dtype=np.int32)
        else:
            k2[ievent]=int(k2[ievent])
            if k2[ievent]<=0 or k2[ievent] > hyper.axes[1].n:
                raise Exception("Illegal value of k2=%d must be between 1 and %d"%(k2[ievent],hyper.axes[1].n))
            self._k2=np.asarray([k2[ievent]-1],dtype=np.int32)

        if ievent >= len(k3):
            self._k3=np.asarray(range(hyper.axes[2].n),dtype=np.int32)
        else:
            k3[ievent]=int(k3[ievent])
            if k3[ievent]<=0 or k3[ievent] > hyper.axes[2].n:
                raise Exception("Illegal value of k3=%d must be between 1 and %d"%(k3[ievent],hyper.axes[2].n))
            self._k3=np.asarray([k3[ievent]-1],dtype=np.int32)

        if ievent >= len(k4):
            self._k4=np.asarray(range(hyper.axes[3].n),dtype=np.int32)
        else:
            k4[ievent]=int(k4[ievent])
            if k4[ievent]<=0 or k4[ievent] > hyper.axes[3].n:
                raise Exception("Illegal value of k4=%d must be between 1 and %d"%(k4[ievent],hyper.axes[3].n))
            self._k4=np.asarray([k4[ievent]-1],dtype=np.int32)
        if ievent >= len(k5):
            self._k5=np.asarray(range(hyper.axes[4].n),dtype=np.int32)
        else:
            k5[ievent]=int(k5[ievent])

            if k5[ievent]<=0 or k5[ievent] > hyper.axes[4].n:
                raise Exception("Illegal value of k5=%d must be between 1 and %d"%(k5[ievent],hyper.axes[4].n))
            self._k5=np.asarray([k5[ievent]-1],dtype=np.int32)

        if ievent >= len(k6):
            self._k6=np.asarray(range(hyper.axes[5].n),dtype=np.int32)
        else:
            k6[ievent]=int(k6[ievent])
            if k6[ievent]<=0 or k6[ievent] > hyper.axes[5].n:
                raise Exception("Illegal value of k6=%d must be between 1 and %d"%(k6[ievent],hyper.axes[5].n))
            self._k6=np.asarray([k6[ievent]-1],dtype=np.int32)




        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make delta functions and impulsive plane waves ')
    parser.add_argument('output', metavar='Output', type=str,
                        help='Output file')                   
    parser.add_argument("--io", type=str,choices=[@GEN_IO_TYPES@], help='IO type. Defaults to defaultIO')
    parser.add_argument("--storage",type=str,choices=["dataByte","dataInt","dataFloat","dataComplex","dataShort",
    "dataComplexDouble","dataDouble"],default="dataFloat")
    parser.add_argument("--memory",type=float,help="Memory in terms of GB",default=.2)
    parser.add_argument("--print_pct",type=float,help="Print progress every X pct (above 100 means no printing)",default=101)
    parser.add_argument("--k1",type=int,help="Fortran index for delta function (first axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--k2",type=int,help="Fortran index for delta function (second axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--k3",type=int,help="Fortran index for delta function (third axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--k4",type=int,help="Fortran index for delta function (fourth axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--k5",type=int,help="Fortran index for delta function (fift axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--k6",type=int,help="Fortran index for delta function (six axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--n1",type=int,help="Number of elements first axis",default=1)
    parser.add_argument("--n2",type=int,help="Number of elements second axis",default=1)
    parser.add_argument("--n3",type=int,help="Number of elements third axis",default=1)
    parser.add_argument("--n4",type=int,help="Number of elements fourth axis",default=1)
    parser.add_argument("--n5",type=int,help="Number of elements fifth axis",default=1)
    parser.add_argument("--n6",type=int,help="Number of elements sixth axis",default=1)
    parser.add_argument("--o1",type=float,help="First location first axis",default=0.)
    parser.add_argument("--o2",type=float,help="First location second axis",default=0.)
    parser.add_argument("--o3",type=float,help="First location third axis",default=0.)
    parser.add_argument("--o4",type=float,help="First location fourth axis",default=0.)
    parser.add_argument("--o5",type=float,help="First location fifth axis",default=0.)
    parser.add_argument("--o6",type=float,help="First location sixth axis",default=0.)
    parser.add_argument("--d1",type=float,help="Sampling first axis",default=1.)
    parser.add_argument("--d2",type=float,help="Sampling second axis",default=1.)
    parser.add_argument("--d3",type=float,help="Sampling third axis",default=1.)
    parser.add_argument("--d4",type=float,help="Sampling fourth axis",default=1.)
    parser.add_argument("--d5",type=float,help="Sampling fifth axis",default=1.)
    parser.add_argument("--d6",type=float,help="Sampling sixth axis",default=1.)
    parser.add_argument("--label1",type=str,help="Label first axis",default="")
    parser.add_argument("--label2",type=str,help="Label second axis",default="")
    parser.add_argument("--label3",type=str,help="Label third axis",default="")
    parser.add_argument("--label4",type=str,help="Label fourth axis",default="")
    parser.add_argument("--label5",type=str,help="Label fifth axis",default="")
    parser.add_argument("--label6",type=str,help="Label sixth axis",default="")
    parser.add_argument("--magnitude",required=True,type=list,help="Magnitude of spikes")

    args = parser.parse_args()
    ioOut=genericIO.defaultIO

    if args.io:
        ioOut=genericIO.io(args.io)

    hyper=hyperFromArgs(args)
    events=[]
    for i in range(len(args.magnitude)):
        events.append(event(i,hyper,args.storage,args.k1,args.k2,args.k3,args.k4,args.k5,args.k6,args.magnitude))

    outFile=genericIO.regFile(ioOut,args.output,storage=args.storage,fromHyper=hyper)
    job=spikeJob(outFile.getStorageType(),events)
    job.setOutputFile(outFile)
    job.setCompleteHyperOut(outFile.getHyper())
    split=GenSplit.serialRegSpace(job, args.memory)
    split.loop(args.print_pct)









